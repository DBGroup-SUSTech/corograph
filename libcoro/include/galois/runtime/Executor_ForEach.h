#ifndef COROGRAPH_EXECUTOR_FOREACH_H
#define COROGRAPH_EXECUTOR_FOREACH_H

#include <algorithm>
#include <functional>
#include <memory>
#include <utility>

#include "galois/runtime/Context.h"
#include "galois/runtime/Corobj.h"

namespace galois {
namespace runtime {

template <class WorkListTy, class Func, class Indexer, class Graph>
class PriorityExecutor {
protected:
  typedef typename WorkListTy::frt_type frt_type;
  typedef typename WorkListTy::msg_type msg_type;
  struct ThreadLocalData {
    std::deque<frt_type> frtbuf;
    std::vector<std::vector<msg_type>> msgbuf;
    std::vector<galois::optional<msg_type>> low, high;
    // the frontier pass the filter.
    std::vector<galois::optional<frt_type>> frt;

    Corobj<bool> coro_sca, coro_gat, coro_syn;
    explicit ThreadLocalData(uint32 numb) : frtbuf(), msgbuf(numb) {
      for (uint32 i = 0; i < numb; i++) {
        // seem not useful, the bottle neck is still memory bound
        msgbuf[i].reserve(100000);
      }
    }
  };

  substrate::TerminationDetection &term;
  WorkListTy wl;
  Func func;
  Graph &graph;

  void doScatter(auto *&frt, auto &tld, bool &didWork) {
    while (frt) {
      didWork = true;
      galois::optional<frt_type> item;
      for (uint32 i = 0; i < frt->count; i++) { // filter
        item = frt->extract_at(i);
        if (func.filterFunc((*item).vid, (*item).dist)) {
          continue;
        }
        // wl.update_cnt("edgecnt", graph.deg[(*item).vid]);
        tld.frt.emplace_back(item);
      }
      frt->clear();
      while (!tld.coro_sca()) {
      }
      tld.coro_syn();
      tld.frt.clear();
      frt = wl.tryPop();
    }
  }

  void doSync(auto &tld) {
    wl.scatter(tld.msgbuf);
    wl.sync();
  }

  void doGather(auto &tld) {
    while (wl.tryPopBlock()) {
      auto *msg = wl.popMsg();
      galois::optional<msg_type> it;
      while (msg) {
        for (uint32 i = 0; i < msg->count; i++) {
          it = msg->extract_at(i);
          if (!((*it).e >> MAX_OFS))
            tld.low.emplace_back(it);
          else
            tld.high.emplace_back(it);
        }
        msg->clear();
        while (!tld.coro_gat()) {
        }
        tld.low.clear();
        tld.high.clear();
        msg = wl.popMsg();
      }
      if (!tld.frtbuf.empty()) {
        wl.push(tld.frtbuf.begin(), tld.frtbuf.end());
        tld.frtbuf.clear();
      }
    }
  }

  bool execute(ThreadLocalData &tld) { // perthread execute
    bool didWork = false;
    auto *f = wl.pop(); // pop a frontier
    while (f) {
      auto start = std::chrono::high_resolution_clock::now();
      doScatter(f, tld, didWork);
      auto end1 = std::chrono::high_resolution_clock::now();
      doSync(tld);
      auto end2 = std::chrono::high_resolution_clock::now();
      doGather(tld);
      auto end3 = std::chrono::high_resolution_clock::now();
      f = wl.slowPop(); // pop frontier from the next priority queue
      auto end4 = std::chrono::high_resolution_clock::now();
      wl.updateTime("1Scatter", start, end1);
      wl.updateTime("2Sync", end1, end2);
      wl.updateTime("3Gather", end2, end3);
      wl.updateTime("4Pop", end3, end4);
    }
    return didWork;
  }

  template <bool isLeader> void go() {
    ThreadLocalData tld(graph.pnum);
    tld.coro_gat = coroGather(tld.frtbuf, tld.low, tld.high);
    tld.coro_sca = coroScatter(tld.msgbuf, tld.frt);
    tld.coro_syn = coroSync(tld.msgbuf);
    do {
      bool b = execute(tld);
      term.localTermination(b);
    } while (!term.globalTermination());
  }

  Corobj<bool> coroScatter(auto &msgbuf, auto &frt) {
    co_yield false;
    for (;;) {
      for (size_t id = 0; id < frt.size(); id += 64) {
        for (size_t prid = id; prid < std::min(id + 64, frt.size()); prid++) {
          _mm_prefetch(&graph.plgraph[(*frt[prid]).vid], _MM_HINT_T0);
        }
        co_yield false;
        for (size_t prid = id; prid < std::min(id + 64, frt.size()); prid++) {
          // scatter
          auto &vtxA = graph.plgraph[(*frt[prid]).vid];
          auto *Arr = vtxA.PE;
          if (vtxA.deg2)
            _mm_prefetch(&graph.pledge[vtxA.offset], _MM_HINT_T0);
          for (uint32 e = 0; e < vtxA.deg1; e += 2) {
            if (!(Arr[e] >> MAX_OFS))
              msgbuf[Arr[e] >> 18].emplace_back(Arr[e], Arr[e + 1],
                                                (*frt[prid]).dist);
            else
              msgbuf[Arr[e] & DW_HALF].emplace_back(
                  MAX_BIT | ((Arr[e] & MIN_HFUP) >> HF_OFST), Arr[e + 1],
                  (*frt[prid]).dist);
          }
          Arr = graph.pledge;
          for (uint32 e = vtxA.offset; e < vtxA.offset + vtxA.deg2; e += 2) {
            if (!(Arr[e] >> MAX_OFS))
              msgbuf[Arr[e] >> 18].emplace_back(Arr[e], Arr[e + 1],
                                                (*frt[prid]).dist);
            else
              msgbuf[Arr[e] & DW_HALF].emplace_back(
                  MAX_BIT | ((Arr[e] & MIN_HFUP) >> HF_OFST), Arr[e + 1],
                  (*frt[prid]).dist);
          }
        }
      }
      co_yield true;
    }
  }
  Corobj<bool> coroGather(auto &frtbuf, auto &low, auto &high) {
    co_yield false;
    for (;;) {
      for (size_t id = 0; id < low.size(); id += 64) {
        for (size_t prid = id; prid < std::min(id + 64, low.size()); prid++) {
          _mm_prefetch(&func.vdata[(*low[prid]).e], _MM_HINT_T0);
        }
        co_yield false;
        for (size_t prid = id; prid < std::min(id + 64, low.size()); prid++) {
          uint32 dst = (*low[prid]).e;
          auto update = func.applyWeight((*low[prid]).w, (*low[prid]).val);
          if (func.gatherFunc(update, dst)) {
            frtbuf.push_back(func.pushFunc(dst, update));
          }
        }
      }
      for (size_t id = 0; id < high.size(); id += 64) {
        for (size_t prid = id; prid < std::min(id + 64, high.size()); prid++) {
          _mm_prefetch(&graph.highedge[(*high[prid]).w], _MM_HINT_T0);
        }
        co_yield false;
        for (size_t prid = id; prid < std::min(id + 64, high.size()); prid++) {
          uint32 deg = (*high[prid]).e & DW_HALF;
          auto dis = (*high[prid]).val;
          for (uint32 e = (*high[prid]).w; e < (*high[prid]).w + deg; e += 2) {
            uint32 dst = graph.highedge[e]; // here
            auto update = func.applyWeight(graph.highedge[e + 1], dis);
            if (func.gatherFunc(update, dst)) {
              frtbuf.push_back(func.pushFunc(dst, update));
            }
          }
        }
      }
      co_yield true;
    }
  }
  Corobj<bool> coroSync(auto &msgbuf) {
    co_yield false;
    for (;;) {
      wl.scatter(msgbuf);
      co_yield true;
    }
  }

public:
  PriorityExecutor(Func _func, Indexer indexer, Graph &_graph)
      : term(substrate::getSystemTermination(activeThreads)),
        wl(_graph.pnum, indexer), func(_func), graph(_graph) {}

  template <typename Frontier> void initFrontier(Frontier &frt) {
    for (auto &f : frt) {
      wl.push(f);
    }
    frt.clear();
  }

  void print_profiling() {
    wl.summary_cnt();
    wl.summaryTime();
  }

  template <typename RangeTy> void initThread(const RangeTy &range) {
    auto rp = range.local_pair();
    wl.push(rp.first, rp.second);
    term.initializeThread();
  }

  void initThread() { term.initializeThread(); }

  void operator()() {
    bool isLeader = substrate::ThreadPool::isLeader();
    if (isLeader)
      go<true>();
    else
      go<false>();
  }
};

template <typename OBIM, typename Graph, typename Indexer, typename Func,
          typename Range> // RangeTy = LocalRange<InsertBag>
void asyncPriorityEdgeMap(Graph &graph, const Indexer &indexer, Func func,
                          const Range &range) {
  auto &barrier = getBarrier(activeThreads);
  PriorityExecutor<OBIM, Func, Indexer, Graph> W(func, indexer, graph);
  substrate::getThreadPool().run(
      activeThreads, [&W, &range]() { W.initThread(range()); },
      std::ref(barrier), std::ref(W));
  W.print_profiling();
}

template <class WorkListTy, class Func, class Indexer, class Graph>
class Executor {
protected:
  typedef typename WorkListTy::frt_type frt_type;
  typedef typename WorkListTy::msg_type msg_type;
  struct ThreadLocalData {
    std::deque<frt_type> frtbuf;
    std::vector<std::vector<msg_type>> msgbuf;
    std::vector<galois::optional<msg_type>> low, high;
    std::vector<galois::optional<frt_type>> frt;

    Corobj<bool> coro_sca, coro_gat, coro_syn;
    explicit ThreadLocalData(uint32 numb) : frtbuf(), msgbuf(numb) {
      for (uint32 i = 0; i < numb; i++) {
        // seem not useful, the bottle neck is still memory bound
        msgbuf[i].reserve(100000);
      }
    }
  };

  substrate::TerminationDetection &term;
  WorkListTy wl;
  Func func;
  Graph &graph;

  void doScatter(auto *&frt, auto &tld, bool &didWork) {
    while (frt) {
      didWork = true;
      galois::optional<frt_type> item;
      for (uint32 i = 0; i < frt->count; i++) { // filter
        item = frt->extract_at(i);
        if (func.filterFunc((*item).vid, (*item).dist)) {
          continue;
        }
        tld.frt.emplace_back(item);
      }
      frt->clear();
      while (!tld.coro_sca())
        tld.coro_syn();
      tld.frt.clear();
      frt = wl.pop();
    }
  }

  void doSync(auto &tld) {
    wl.scatter(tld.msgbuf);
    wl.sync();
  }

  void doGather(auto &tld) {
    while (wl.tryPopBlock()) {
      auto *msg = wl.popMsg();
      galois::optional<msg_type> it;
      while (msg) {
        for (uint32 i = 0; i < msg->count; i++) {
          it = msg->extract_at(i);
          if (!((*it).e >> MAX_OFS))
            tld.low.emplace_back(it);
          else
            tld.high.emplace_back(it);
        }
        msg->clear();
        while (!tld.coro_gat()) {
          if (!tld.frtbuf.empty()) {
            wl.push(tld.frtbuf.begin(), tld.frtbuf.end());
            tld.frtbuf.clear();
          }
        }
        tld.low.clear();
        tld.high.clear();
        msg = wl.popMsg();
      }
    }
  }

  bool execute(ThreadLocalData &tld) { // perthread execute
    bool didWork = false;
    auto *p = wl.pop(); // pop a chunk
    while (p) {
      auto start = std::chrono::high_resolution_clock::now();
      doScatter(p, tld, didWork);
      auto end1 = std::chrono::high_resolution_clock::now();
      doSync(tld);
      auto end2 = std::chrono::high_resolution_clock::now();
      doGather(tld);
      auto end3 = std::chrono::high_resolution_clock::now();
      p = wl.slowPop(); // pop the next priority queue
      wl.updateTime("1Scatter", start, end1);
      wl.updateTime("2Sync", end1, end2);
      wl.updateTime("3Gather", end2, end3);
    }
    return didWork;
  }

  template <bool isLeader> void go() {
    ThreadLocalData tld(graph.pnum);
    tld.coro_gat = coroGather(tld.frtbuf, tld.low, tld.high);
    tld.coro_sca = coroScatter(tld.msgbuf, tld.frt);
    tld.coro_syn = coroSync(tld.msgbuf);
    do {
      bool b = execute(tld);
      term.localTermination(b);
    } while (!term.globalTermination());
  }

  Corobj<bool> coroScatter(auto &msgbuf, auto &frt) {
    co_yield false;
    for (;;) {
      for (size_t id = 0; id < frt.size(); id += 64) {
        for (size_t prid = id; prid < std::min(id + 64, frt.size()); prid++) {
          _mm_prefetch(&graph.plgraph[(*frt[prid]).vid], _MM_HINT_T0);
        }
        co_yield false;
        for (size_t prid = id; prid < std::min(id + 64, frt.size()); prid++) {
          // scatter
          auto &vtxA = graph.plgraph[(*frt[prid]).vid];
          auto *Arr = vtxA.PE;
          if (vtxA.deg2)
            _mm_prefetch(&graph.pledge[vtxA.offset], _MM_HINT_T0);
          for (uint32 e = 0; e < vtxA.deg1; e += 2) {
            if (!(Arr[e] >> MAX_OFS))
              msgbuf[Arr[e] >> 18].emplace_back(Arr[e], Arr[e + 1],
                                                (*frt[prid]).dist);
            else
              msgbuf[Arr[e] & DW_HALF].emplace_back(
                  MAX_BIT | ((Arr[e] & MIN_HFUP) >> HF_OFST), Arr[e + 1],
                  (*frt[prid]).dist);
          }
          Arr = graph.pledge;
          for (uint32 e = vtxA.offset; e < vtxA.offset + vtxA.deg2; e += 2) {
            if (!(Arr[e] >> MAX_OFS))
              msgbuf[Arr[e] >> 18].emplace_back(Arr[e], Arr[e + 1],
                                                (*frt[prid]).dist);
            else
              msgbuf[Arr[e] & DW_HALF].emplace_back(
                  MAX_BIT | ((Arr[e] & MIN_HFUP) >> HF_OFST), Arr[e + 1],
                  (*frt[prid]).dist);
          }
        }
      }
      co_yield true;
    }
  }
  Corobj<bool> coroGather(auto &frtbuf, auto &low, auto &high) {
    co_yield false;
    for (;;) {
      for (size_t id = 0; id < low.size(); id += 64) {
        for (size_t prid = id; prid < std::min(id + 64, low.size()); prid++) {
          _mm_prefetch(&func.vdata[(*low[prid]).e], _MM_HINT_T0);
        }
        co_yield false;
        for (size_t prid = id; prid < std::min(id + 64, low.size()); prid++) {
          uint32 dst = (*low[prid]).e;
          auto update = func.applyWeight((*low[prid]).w, (*low[prid]).val);
          if (func.gatherFunc(update, dst)) {
            frtbuf.push_back(func.pushFunc(dst, update));
          }
        }
      }
      for (size_t id = 0; id < high.size(); id += 64) {
        for (size_t prid = id; prid < std::min(id + 64, high.size()); prid++) {
          _mm_prefetch(&graph.highedge[(*high[prid]).w], _MM_HINT_T0);
        }
        co_yield false;
        for (size_t prid = id; prid < std::min(id + 64, high.size()); prid++) {
          uint32 deg = (*high[prid]).e & DW_HALF;
          auto dis = (*high[prid]).val;
          for (uint32 e = (*high[prid]).w; e < (*high[prid]).w + deg; e += 2) {
            uint32 dst = graph.highedge[e]; // here
            auto update = func.applyWeight(graph.highedge[e + 1], dis);
            if (func.gatherFunc(update, dst)) {
              frtbuf.push_back(func.pushFunc(dst, update));
            }
          }
        }
      }
      co_yield true;
    }
  }
  Corobj<bool> coroSync(auto &msgbuf) {
    co_yield false;
    for (;;) {
      wl.scatter(msgbuf);
      co_yield true;
    }
  }

public:
  Executor(Func _func, Indexer indexer, Graph &_graph)
      : term(substrate::getSystemTermination(activeThreads)),
        wl(_graph.pnum, indexer), func(_func), graph(_graph) {}

  template <typename Frontier> void initFrontier(Frontier &frt) {
    for (auto &a : frt) {
      wl.push(a);
    }
    frt.clear();
  }

  void print_profiling() {
    wl.summary_cnt();
    wl.summaryTime();
  }

  void initThread() { term.initializeThread(); }

  void operator()() {
    bool isLeader = substrate::ThreadPool::isLeader();
    if (isLeader)
      go<true>();
    else
      go<false>();
  }
};
template <typename OBIM, typename Graph, typename Frontier, typename Indexer,
          typename Func> // RangeTy = LocalRange<InsertBag>
void syncEdgeMap(Graph &graph, Frontier &frontier, const Indexer &indexer,
                 Func func) {
  auto &barrier = getBarrier(activeThreads);
  Executor<OBIM, Func, Indexer, Graph> W(func, indexer, graph);
  W.initFrontier(frontier);
  substrate::getThreadPool().run(activeThreads, [&W]() { W.initThread(); });
  substrate::getThreadPool().run(activeThreads, std::ref(W));
  W.print_profiling();
}
} // end namespace runtime
} // end namespace galois
#endif
