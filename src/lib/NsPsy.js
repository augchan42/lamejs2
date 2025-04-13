import * as common from './common.js';
const { System, VbrMode, Float, ShortBlock, Util, Arrays, new_array_n, new_byte, new_double, new_float, new_float_n, new_int, new_int_n, assert } = common;

import { Encoder } from './Encoder.js';

/**
 * Variables used for --nspsytune
 *
 * @author Ken
 */
class NsPsy {
    constructor() {
        this.last_en_subshort = new_float_n([4, 9]);
        this.lastAttacks = new_int(4);
        this.pefirbuf = new_float(19);
        this.longfact = new_float(Encoder.SBMAX_l);
        this.shortfact = new_float(Encoder.SBMAX_s);

        /**
         * short block tuning
         */
        this.attackthre = 0.;
        this.attackthre_s = 0.;
    }
}

export { NsPsy };
