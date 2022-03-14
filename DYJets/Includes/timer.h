#ifndef TIMER_H
#define TIMER_H

#include <atomic>
#include <thread>

namespace util
{

/**
 * \brief Displays progress in the terminal.
 *
 * The \c timer class displays progress for a process consisting of a given
 * amount of steps. It will keep updating even if the calling program is
 * blocked.
 *
 * \note Do to the multi-threaded nature of this class, you may experience
 *       interleaved output. There is no way around it besides explicit locking,
 *       whose implementation is left as an exercise to the reader.
 *
 * \note This class isn't thread-safe nor reentrant.
 */
class timer
{
    using clock = std::chrono::steady_clock;
    using time_point = clock::time_point;
    using duration = clock::duration;

  public:
    using counter = unsigned long;

  private:
    // In both threads
    const counter _steps, _from;
    std::atomic<counter> _done;
    std::atomic<bool> _stop;

    // In control thread
    bool _running = false;
    std::thread _thread;

    // In timer thread
    time_point _start, _last;
    counter _lastcounter = 0;

  public:
    /**
     * \brief Creates a new timer that will count up to the given \c last step,
     *        starting from \c begin.
     */
    explicit timer(counter last, counter begin = 0);

    /** \brief Destructor. */
    virtual ~timer();

    /**
     * \brief Starts counting.
     */
    void start();

    /**
     * \brief Stops counting.
     * \note Resuming counting after it is stopped isn't supported.
     */
    void stop();

    /** \brief Returns the total number of steps to be performed. */
    counter steps() const { return _steps; };

    /** \brief Returns the number of steps already done. */
    counter done() const { return _done; };

    /**
     * \brief Sets the number of steps already done.
     * \see next()
     */
    void update(counter done) { _done = done; };

    /**
     * \brief Increases the number of steps already done by one.
     * \see update()
     */
    void next() { ++_done; };

    /**
     * \brief Check whether the timer is running.
     */
    bool isrunning() const;

  private:
    void displayprogress(const time_point &now, bool at_end = false) const;
};
} // namespace util

#endif // TIMER_H
