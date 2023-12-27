#include <atomic>
#include <coroutine>
#include <thread>

inline static void asmPause() { asm volatile("pause"); }

class CoroLock {
public:
  CoroLock() : flag(false) {}
  void lock(std::coroutine_handle<> hd) {
    while (flag.test_and_set(std::memory_order_acquire)) {
      if (!hd.done())
        hd();
      else
        std::this_thread::yield();
    }
  }
  void lock() {
    while (flag.test_and_set(std::memory_order_acquire)) {
      asmPause();
    }
  }
  void unlock() { flag.clear(std::memory_order_release); }

private:
  std::atomic_flag flag;
};
