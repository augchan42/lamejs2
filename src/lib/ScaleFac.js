//package mp3;

/**
 * Layer III side information.
 *
 * @author Ken
 *
 */

import * as common from './common.js';
const { System, VbrMode, Float, ShortBlock, Util, Arrays, new_array_n, new_byte, new_double, new_float, new_float_n, new_int, new_int_n, assert } = common;

import { Encoder } from './Encoder.js';

class ScaleFac {
    constructor(arrL, arrS, arr21, arr12) {
        this.l = new_int(1 + Encoder.SBMAX_l);
        this.s = new_int(1 + Encoder.SBMAX_s);
        this.psfb21 = new_int(1 + Encoder.PSFB21);
        this.psfb12 = new_int(1 + Encoder.PSFB12);
        const l = this.l;
        const s = this.s;

        if (arguments.length == 4) {
            this.arrL = arrL;
            this.arrS = arrS;
            this.arr21 = arr21;
            this.arr12 = arr12;

            System.arraycopy(this.arrL, 0, l, 0, Math.min(this.arrL.length, this.l.length));
            System.arraycopy(this.arrS, 0, s, 0, Math.min(this.arrS.length, this.s.length));
            System.arraycopy(this.arr21, 0, this.psfb21, 0, Math.min(this.arr21.length, this.psfb21.length));
            System.arraycopy(this.arr12, 0, this.psfb12, 0, Math.min(this.arr12.length, this.psfb12.length));
        }
    }
}

export { ScaleFac };
