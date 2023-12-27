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

#ifndef GALOIS_WORKLIST_OBIM_H
#define GALOIS_WORKLIST_OBIM_H

#include <deque>
#include <limits>
#include <type_traits>

#include "galois/FlatMap.h"
#include "galois/runtime/Substrate.h"
#include "galois/substrate/PerThreadStorage.h"
#include "galois/substrate/Termination.h"
#include "galois/worklists/Chunk.h"
#include "galois/worklists/WorkListHelpers.h"
#include "utils.h"

namespace galois {
namespace worklists {

template <typename Indexer, typename PQueue, typename Block, typename Frontier,
          typename Message>
struct OBIM {

  typedef typename PQueue::frt_type frt_type;
  typedef typename Indexer::compare Compare;
  typedef typename Block::msg_type msg_type;
  typedef uint32 Index;

  Compare compare;
  Index earliest;

  struct ThreadData {
    galois::flat_map<Index, PQueue *, std::less<Index>> local;
    Index curIndex;
    Index scanStart;
    PQueue *current;
    unsigned int lastMasterVersion;
    unsigned int numB;
    Block *currentB; // current owned block pointer for gather
    uint32 Bid;
    Block **localBP; // local block pointer, for thread local

    // for profiling only
    std::map<std::string, std::map<uint32, uint32>> cnt;
    std::map<std::string, uint64_t> cntTime;

    ThreadData(Index initial, uint32 _numB)
        : curIndex(initial), scanStart(initial), current(0), currentB(0),
          lastMasterVersion(0), numB(_numB) {
      localBP = new Block *[_numB];
    }

    void updateBP(Block **BP) {
      for (uint32 i = 0; i < numB; i++) {
        localBP[i] = BP[i];
      }
    }
  };

  substrate::PerThreadStorage<ThreadData> data;
  substrate::PaddedLock<true> masterLock;
  std::deque<std::pair<Index, PQueue *>> masterLog;

  // The list of block ids need to be gathered.
  substrate::PerSocketStorage<ConExtLinkedQueue<Block, true>> blockList;
  Block **blockPointer; // per block pointer
  uint32 numB;          // block number

  std::atomic<unsigned int> masterVersion;
  Indexer indexer;

  bool updateLocal(ThreadData &p) {
    if (p.lastMasterVersion != masterVersion.load(std::memory_order_relaxed)) {
      for (;
           p.lastMasterVersion < masterVersion.load(std::memory_order_relaxed);
           ++p.lastMasterVersion) {
        std::pair<Index, PQueue *> logEntry = masterLog[p.lastMasterVersion];
        p.local[logEntry.first] = logEntry.second;
      }
      return true;
    }
    return false;
  }

  __attribute__((noinline)) PQueue *slowUpdateLocalOrCreate(ThreadData &p,
                                                            Index i) {
    do {
      updateLocal(p);
      auto it = p.local.find(i);
      if (it != p.local.end())
        return it->second;
    } while (!masterLock.try_lock());
    updateLocal(p);
    auto it = p.local.find(i);
    PQueue *pq = (it != p.local.end()) ? it->second : nullptr;
    if (!pq) {
      pq = new PQueue();
      p.local[i] = pq;
      p.lastMasterVersion = masterVersion.load(std::memory_order_relaxed) + 1;
      masterLog.push_back(std::make_pair(i, pq));
      masterVersion.fetch_add(1);
    }
    masterLock.unlock();
    return pq;
  }

  inline PQueue *updateLocalOrCreate(ThreadData &p, Index i) {
    auto it = p.local.find(i);
    if (it != p.local.end())
      return it->second;
    return slowUpdateLocalOrCreate(p, i);
  }

public:
  OBIM(uint32 _numB, const Indexer &x = Indexer())
      : numB(_numB), masterVersion(0), indexer(x), data(this->earliest, _numB),
        earliest(std::numeric_limits<Index>::min()) {
    blockPointer = new Block *[_numB];
    for (uint32 i = 0; i < _numB; i++) {
      blockPointer[i] = new Block(i);
    }
    for (unsigned i = 0; i < runtime::activeThreads; ++i) {
      ThreadData &o = *data.getRemote(i);
      o.updateBP(blockPointer);
    }
  }

  ~OBIM() {
    for (auto ii = masterLog.rbegin(), ei = masterLog.rend(); ii != ei; ++ii) {
      delete ii->second;
    }
  }

  uint32 getBid() {
    ThreadData &p = *data.getLocal();
    p.Bid = p.currentB->Bid;
    return p.Bid;
  }

  void updateTime(std::string name, auto &st, auto &ed) {
    auto duration =
        std::chrono::duration_cast<std::chrono::nanoseconds>(ed - st);
    ThreadData &p = *data.getLocal();
    p.cntTime[name] += duration.count();
  }
  void summaryTime() {
    uint32 thd = runtime::activeThreads;
    std::map<std::string, std::vector<uint64>> sumTime;
    for (unsigned i = 0; i < runtime::activeThreads; ++i) {
      auto &pcnt = data.getRemote(i)->cntTime;
      for (auto &name : pcnt)
        sumTime[name.first].push_back(name.second);
    }
    uint64_t totaltime = 0;
    for (auto &name : sumTime) {
      uint64_t tsum = 0;
      for (auto t : sumTime[name.first]) {
        tsum += t;
        printf("%lu ", t / 1000000);
      }
      totaltime += tsum;
      printf("\n  %s sumtime: %lu ms\n", name.first.c_str(), tsum / 1000000);
    }
    printf("totalTime: %lu avgTime: %lu\n", totaltime / 1000000,
           totaltime / 1000000 / getActiveThreads());
  }

  void update_cnt(std::string name, uint32 val, uint32 num = 1) {
    ThreadData &p = *data.getLocal();
    p.cnt[name][val] += num;
  }

  void summary_cnt() {
    std::map<std::string, std::map<uint32, uint32>> sum;
    for (unsigned i = 0; i < runtime::activeThreads; ++i) {
      auto &pcnt = data.getRemote(i)->cnt;
      for (auto &name : pcnt)
        for (auto &val : name.second)
          sum[name.first][val.first] += val.second;
    }
    // uint32 level[] = {1000,65536>>1, 65536, 262144>>1,
    // 262144, 262144<<1, 262144<<2, 0xffffffff};
    uint32 level[] = {10, 100, 1000, 10000, 100000, 1000000};
    uint32 levelnum = 6;
    for (auto &name : sum) {
      uint32 cnt[levelnum], scnt[levelnum];
      for (uint32 i = 0; i < levelnum; i++)
        cnt[i] = scnt[i] = 0;
      printf("summary of %s\n", name.first.c_str());
      for (auto &it : name.second) {
        for (uint32 l = 0; l < levelnum; l++) {
          if (it.first < level[l]) {
            cnt[l] += it.second;
            scnt[l] += it.second * it.first;
            break;
          }
        }
      }
      uint32 sumc = 0, sums = 0;
      for (uint32 i = 0; i < levelnum; i++) {
        sumc += cnt[i];
        sums += scnt[i];
      }
      printf("sumc: %u sums:%u\n", sumc, sums);
      for (uint32 l = 0; l < levelnum; l++) {
        printf("%u-%u: cnt %u sum %u %.2f%% %.2f%%\n",
               (l == 0) ? 0 : level[l - 1], level[l], cnt[l], scnt[l],
               (double)100 * cnt[l] / sumc, (double)100 * scnt[l] / sums);
      }
    }
  }

  void scatter(uint32 bid, const msg_type &msg) {
    ThreadData &p = *data.getLocal();
    scatter(p, bid, msg);
  }
  void scatter(ThreadData &p, uint32 bid, const msg_type &msg) {
    if (p.localBP[bid]->push(msg)) { // empty and not in list
      blockList.getLocal()->push(p.localBP[bid]);
    }
  }

  template <typename Iter>
  void scatter(ThreadData &p, uint32 bid, Iter b, Iter e) {
    while (b != e)
      scatter(p, bid, *b++);
  }

  void scatter(auto &msgbuf) {
    ThreadData &p = *data.getLocal();
    for (uint32 bid = 0; bid < msgbuf.size(); bid++) {
      // update_cnt("msgbufsize", msgbuf[i].size());
      scatter(p, bid, msgbuf[bid].begin(), msgbuf[bid].end());
      // for (auto iter = msgbuf[bid].begin(); iter != msgbuf[bid].end();) {
      //   if (p.localBP[bid]->push(*iter++)) { // empty and not in list
      //     blockList.getLocal()->push(p.localBP[bid]);
      //   }
      // }
      msgbuf[bid].clear();
    }
  }

  void sync() {
    ThreadData &p = *data.getLocal();
    for (uint32 i = 0; i < numB; i++) {
      if (p.localBP[i]->pushNext()) {
        // push all local update to block buffer
        blockList.getLocal()->push(p.localBP[i]);
      }
    }
  }

  Block *popBlockByID(unsigned int i) {
    auto &I = *blockList.getRemote(i);
    return I.pop();
  }

  Block *popBlock() {
    int id = substrate::ThreadPool::getTID(); // thread id
    Block *r = popBlockByID(id);
    if (r)
      return r;

    for (int i = id + 1; i < (int)blockList.size(); ++i) {
      r = popBlockByID(i);
      if (r)
        return r;
    }
    for (int i = 0; i < id; ++i) {
      r = popBlockByID(i);
      if (r)
        return r;
    }
    return 0;
  }

  bool tryPopBlock() { // try to pop a block from blockList
    ThreadData &p = *data.getLocal();
    p.currentB = popBlock(); // blockList.getLocal()->pop();
    if (p.currentB) {
      return true;
    }
    return false;
  }

  Message *popMsg() {
    ThreadData &p = *data.getLocal();
    Message *item = p.currentB->pop();
    return item;
  }

  void push(const frt_type &val) {
    Index index = indexer(val);
    ThreadData &p = *data.getLocal();

    if (index == p.curIndex && p.current) {
      p.current->push(val);
      return;
    }

    PQueue *pq = updateLocalOrCreate(p, index);
    if (this->compare(index, p.scanStart))
      p.scanStart = index;
    if (this->compare(index, p.curIndex)) {
      p.curIndex = index;
      p.current = pq;
    }
    pq->push(val);
  }

  template <typename Iter> void push(Iter b, Iter e) {
    while (b != e)
      push(*b++);
  }

  Frontier *pop() { // pop a frontier == tryPop + slowPop
    ThreadData &p = *data.getLocal();
    PQueue *pq = p.current;
    Frontier *item = nullptr;
    if (pq)
      item = pq->pop();
    if (item)
      return item;
    return slowPop();
  }

  Frontier *tryPop() { // try to pop a frontier
    ThreadData &p = *data.getLocal();
    PQueue *pq = p.current;
    Frontier *item = nullptr;
    if (pq)
      item = pq->pop();
    return item;
  }

  __attribute__((noinline)) Frontier *slowPop() {
    ThreadData &p = *data.getLocal();
    bool localLeader = substrate::ThreadPool::isLeader();
    Index msS = this->earliest;
    updateLocal(p);
    msS = p.scanStart;
    if (localLeader) {
      for (unsigned i = 0; i < runtime::activeThreads; ++i) {
        Index o = data.getRemote(i)->scanStart;
        if (this->compare(o, msS))
          msS = o;
      }
    } else {
      Index o = data.getRemote(substrate::ThreadPool::getLeader())->scanStart;
      if (this->compare(o, msS))
        msS = o;
    }
    for (auto ii = p.local.lower_bound(msS), ei = p.local.end(); ii != ei;
         ++ii) {
      Frontier *item = ii->second->pop();
      if (item) {
        p.current = ii->second;
        p.curIndex = ii->first;
        p.scanStart = ii->first;
        return item;
      }
    }
    return nullptr;
  }
};
} // end namespace worklists
} // end namespace galois

#endif
