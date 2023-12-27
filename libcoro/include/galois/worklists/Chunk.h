/*
 * This file belongs to the Galois project, a C++ library for exploiting
 * parallelism. The code is being released under the terms of the 3-Clause BSD
 * License (a copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

#ifndef GALOIS_WORKLIST_CHUNK_H
#define GALOIS_WORKLIST_CHUNK_H

#include "galois/FixedSizeRing.h"
#include "galois/config.h"
#include "galois/runtime/Mem.h"
#include "galois/substrate/PaddedLock.h"
#include "galois/worklists/WLCompileCheck.h"
#include "galois/worklists/WorkListHelpers.h"

namespace galois {
namespace runtime {
extern unsigned activeThreads;
}
namespace worklists {

namespace internal {

template <typename T, int ChunkSize>
class Chunk : public FixedSizeRing2<T, ChunkSize> {
  Chunk *next;

public:
  Chunk() : next(0) {}
  Chunk *&getNext() { return next; }
  Chunk *const &getNext() const { return next; }
};

template <typename T, int MessageSize>
struct Block : public ConExtLinkedQueue<Block<T, MessageSize>, true>::ListNode {
public:
  using Message = Chunk<T, MessageSize>;
  runtime::FixedSizeAllocator<Message> alloc;
  struct MsgBuf {
    Message *cur;
    Message *next;
    MsgBuf() : cur(0), next(0) {}
  };
  substrate::PerThreadStorage<MsgBuf> data;
  MessageQueue<Message, true> buffer;
  uint32_t bid;

  Message *mkMessage() {
    Message *ptr = alloc.allocate(1);
    alloc.construct(ptr);
    return ptr;
  }

  void delMessage(Message *ptr) {
    alloc.destroy(ptr);
    alloc.deallocate(ptr, 1);
  }
  bool pushMessage(Message *m) { return buffer.push(m); }
  Message *popMessage() { return buffer.pop(); }

  template <typename... Args> bool emplacei(MsgBuf &msg, Args &&...args) {
    if (msg.next && !msg.next->emplace_back(std::forward<Args>(args)...))
      return false;
    if (msg.next) { // full
      bool ret = pushMessage(msg.next);
      msg.next = 0;
      return ret;
    }
    msg.next = mkMessage();
    msg.next->emplace_back(std::forward<Args>(args)...);
    return false;
  }

  typedef T msg_type;

  explicit Block(uint32_t _bid) : bid(_bid) {}
  Block(const Block &) = delete;
  Block &operator=(const Block &) = delete;

  bool push(const T &val) {
    MsgBuf &msg = *data.getLocal();
    return emplacei(msg, val);
  }

  bool pushNext() {
    bool ret = false;
    MsgBuf &msg = *data.getLocal();
    if (msg.next) { //&& !msg.next->empty()
      ret = pushMessage(msg.next);
      msg.next = 0;
      // msg.next = mkMessage();
    }
    return ret;
  }

  Message *pop() {
    MsgBuf &msg = *data.getLocal();
    // if (msg.cur && !msg.cur->empty()){
    //     return msg.cur;
    // }
    if (msg.cur)
      delMessage(msg.cur);
    msg.cur = popMessage();
    if (!msg.cur) {
      msg.cur = msg.next;
      msg.next = 0;
    }
    return msg.cur;
    // if (msg.cur && !msg.cur->empty())
    //     return msg.cur;
    // return nullptr;
  }
};

template <typename T, int FrontierSize> struct PQueue {
public:
  using Frontier = Chunk<T, FrontierSize>;
  runtime::FixedSizeAllocator<Frontier> alloc;

  struct FrtBuf {
    Frontier *cur;
    Frontier *next;
    FrtBuf() : cur(0), next(0) {}
  };

  substrate::PerThreadStorage<FrtBuf> data;
  substrate::PerSocketStorage<ConExtLinkedQueue<Frontier, true>> Q;

  Frontier *mkFrontier() {
    Frontier *ptr = alloc.allocate(1);
    alloc.construct(ptr);
    return ptr;
  }

  void delFrontier(Frontier *ptr) {
    alloc.destroy(ptr);
    alloc.deallocate(ptr, 1);
  }

  void pushFrontier(Frontier *f) {
    auto &I = *Q.getLocal();
    I.push(f);
  }

  Frontier *popFrontierByID(unsigned int i) {
    auto &I = *Q.getRemote(i);
    return I.pop();
  }

  Frontier *popFrontier() {
    int id = substrate::ThreadPool::getTID(); // thread id
    Frontier *r = popFrontierByID(id);
    if (r)
      return r;

    for (int i = id + 1; i < (int)Q.size(); ++i) {
      r = popFrontierByID(i);
      if (r)
        return r;
    }
    for (int i = 0; i < id; ++i) {
      r = popFrontierByID(i);
      if (r)
        return r;
    }
    return 0;
  }

  template <typename... Args> void emplacei(FrtBuf &frt, Args &&...args) {
    if (frt.next && !frt.next->emplace_back(std::forward<Args>(args)...))
      return;
    if (frt.next) {
      pushFrontier(frt.next);
      frt.next = 0;
      return;
    }
    frt.next = mkFrontier();
    frt.next->emplace_back(std::forward<Args>(args)...);
  }

  typedef T frt_type;

  PQueue() = default;
  PQueue(const PQueue &) = delete;
  PQueue &operator=(const PQueue &) = delete;

  void push(const T &val) {
    FrtBuf &frt = *data.getLocal();
    emplacei(frt, val);
  }

  template <typename Iter> void push(Iter b, Iter e) {
    FrtBuf &frt = *data.getLocal();
    while (b != e)
      emplacei(frt, *b++);
  }

  Frontier *pop() {
    FrtBuf &frt = *data.getLocal();
    // if (frt.cur && !frt.cur->empty()){
    //     return frt.cur;
    // }
    if (frt.cur)
      delFrontier(frt.cur);
    frt.cur = popFrontier();
    if (!frt.cur) {       // if queue is empty
      frt.cur = frt.next; // use next as current
      frt.next = 0;
    }
    return frt.cur;
    // if (frt.cur && !frt.cur->empty())
    //     return frt.cur;
    // return nullptr;
  }

  // galois::optional<value_type> pop() {
  //     p& frt = *data.getLocal();
  //     galois::optional<value_type> retval;
  //     if (frt.cur && (retval = frt.cur->extract_front()))
  //         return retval;
  //     if (frt.cur)
  //         delFrontier(frt.cur);
  //     frt.cur = popFrontier();
  //     if (!frt.cur) {
  //         frt.cur  = frt.next;
  //         frt.next = 0;
  //     }
  //     if (frt.cur)
  //         return frt.cur->extract_front();
  //     return galois::optional<value_type>();
  // }
};

} // namespace internal
template <int ChunkSize = 64, typename T = int>
using PQ = internal::PQueue<T, ChunkSize>;
template <int ChunkSize = 64, typename T = int>
using BL = internal::Block<T, ChunkSize>;
template <int ChunkSize = 64, typename T = int>
using CK = internal::Chunk<T, ChunkSize>;
} // end namespace worklists
} // end namespace galois

#endif
