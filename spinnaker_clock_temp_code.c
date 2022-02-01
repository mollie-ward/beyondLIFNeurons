#include <sark.h>
#include <stdfix-full-iso.h>
#include "debug.h"
#include "spin1_api.h"
#include <stdfix.h>
#include <stdfix-exp.h>
#include "stdfix-exp.h"

uint32_t num_ticks = 10;
uint32_t tick = 0;
uint32_t measure_in, measure_out;

void tick_callback(uint ticks, uint dummy) {
    if (tick < num_ticks) {
        char buf[64];
        measure_in = tc[T1_COUNT];

        // ADD CODE IN HERE THAT NEEDS TO BE TIMED

        
        measure_out = tc[T1_COUNT];
        io_printf(IO_BUF,
            "Iteration: %u    clock in: %u     clock out: %u     total time: %u \n",
            tick, measure_in, measure_out, measure_in - measure_out);
        }
		tick++;
	} else {
		spin1_exit(0);
	}
}

void c_main(void) {
	spin1_set_timer_tick(1000000);  // set timer tick to 1 second (not 1 ms!)
	spin1_callback_on(TIMER_TICK, tick_callback, 1);  // timer callback
	spin1_start(SYNC_NOWAIT);  // start event-driven operation
}
