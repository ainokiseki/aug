#pragma once
#include"metric.h"
#include<vector>
#include<unordered_map>
#include<mutex>
#include<condition_variable>
#include<set>
#include<iostream>

inline int def_ij(int x, int y, int n) {



	return n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;
}

void cacap(Eigen::VectorXd* res, int start, int end, std::vector<int64_t>* aj1, Eigen::VectorXd* rhs,int n);


class ReaderWriterLatch {
    using mutex_t = std::mutex;
    using cond_t = std::condition_variable;


public:
    ReaderWriterLatch() = default;
    ~ReaderWriterLatch() { std::lock_guard<mutex_t> guard(mutex_); }



    /**
     * Acquire a write latch.
     */
    void WLock() {
        std::unique_lock<mutex_t> latch(mutex_);
        while (writer_entered_) {
            reader_.wait(latch);
        }
        writer_entered_ = true;
        while (reader_count_ > 0) {
            writer_.wait(latch);
        }
    }

    /**
     * Release a write latch.
     */
    void WUnlock() {
        std::lock_guard<mutex_t> guard(mutex_);
        writer_entered_ = false;
        reader_.notify_all();
    }

    /**
     * Acquire a read latch.
     */
    void RLock() {
        std::unique_lock<mutex_t> latch(mutex_);
        while (writer_entered_ ) {
            reader_.wait(latch);
        }
        reader_count_++;
    }

    /**
     * Release a read latch.
     */
    void RUnlock() {
        std::lock_guard<mutex_t> guard(mutex_);
        reader_count_--;
        if (writer_entered_) {
            if (reader_count_ == 0) {
                writer_.notify_one();
            }
        }

    }

private:
    mutex_t mutex_;
    cond_t writer_;
    cond_t reader_;
    uint32_t reader_count_{ 0 };
    bool writer_entered_{ false };
};
