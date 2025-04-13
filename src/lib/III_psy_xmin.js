import { Encoder } from './Encoder.js';
import * as common from './common.js';
const { System, VbrMode, Float, ShortBlock, Util, Arrays, new_array_n, new_byte, new_double, new_float, new_float_n, new_int, new_int_n, assert } = common;

class III_psy_xmin {
    constructor() {
        this.l = new_float(Encoder.SBMAX_l);
        this.s = new_float_n([Encoder.SBMAX_s, 3]);
    }

    assign(iii_psy_xmin) {
        System.arraycopy(iii_psy_xmin.l, 0, this.l, 0, Encoder.SBMAX_l);
        for (var i = 0; i < Encoder.SBMAX_s; i++) {
            for (var j = 0; j < 3; j++) {
                this.s[i][j] = iii_psy_xmin.s[i][j];
            }
        }
    }
}

export { III_psy_xmin };
