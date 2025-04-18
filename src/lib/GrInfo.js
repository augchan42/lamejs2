//package mp3;
import * as common from './common.js';
const { System, VbrMode, Float, ShortBlock, Util, Arrays, new_array_n, new_byte, new_double, new_float, new_float_n, new_int, new_int_n, assert } = common;

import { L3Side } from './L3Side.js';

class GrInfo {
    constructor() {
        //float xr[] = new float[576];
        this.xr = new_float(576);
        //int l3_enc[] = new int[576];
        this.l3_enc = new_int(576);
        //int scalefac[] = new int[L3Side.SFBMAX];
        this.scalefac = new_int(L3Side.SFBMAX);
        this.xrpow_max = 0.;

        this.part2_3_length = 0;
        this.big_values = 0;
        this.count1 = 0;
        this.global_gain = 0;
        this.scalefac_compress = 0;
        this.block_type = 0;
        this.mixed_block_flag = 0;
        this.table_select = new_int(3);
        this.subblock_gain = new_int(3 + 1);
        this.region0_count = 0;
        this.region1_count = 0;
        this.preflag = 0;
        this.scalefac_scale = 0;
        this.count1table_select = 0;

        this.part2_length = 0;
        this.sfb_lmax = 0;
        this.sfb_smin = 0;
        this.psy_lmax = 0;
        this.sfbmax = 0;
        this.psymax = 0;
        this.sfbdivide = 0;
        this.width = new_int(L3Side.SFBMAX);
        this.window = new_int(L3Side.SFBMAX);
        this.count1bits = 0;
        /**
         * added for LSF
         */
        this.sfb_partition_table = null;
        this.slen = new_int(4);

        this.max_nonzero_coeff = 0;
    }

    static clone_int(array) {
        return new Int32Array(array);
    }

    static clone_float(array) {
        return new Float32Array(array);
    }

    assign(other) {
        this.xr = GrInfo.clone_float(other.xr);
        this.l3_enc = GrInfo.clone_int(other.l3_enc);
        this.scalefac = GrInfo.clone_int(other.scalefac);
        this.xrpow_max = other.xrpow_max;

        this.part2_3_length = other.part2_3_length;
        this.big_values = other.big_values;
        this.count1 = other.count1;
        this.global_gain = other.global_gain;
        this.scalefac_compress = other.scalefac_compress;
        this.block_type = other.block_type;
        this.mixed_block_flag = other.mixed_block_flag;
        this.table_select = GrInfo.clone_int(other.table_select);
        this.subblock_gain = GrInfo.clone_int(other.subblock_gain);
        this.region0_count = other.region0_count;
        this.region1_count = other.region1_count;
        this.preflag = other.preflag;
        this.scalefac_scale = other.scalefac_scale;
        this.count1table_select = other.count1table_select;

        this.part2_length = other.part2_length;
        this.sfb_lmax = other.sfb_lmax;
        this.sfb_smin = other.sfb_smin;
        this.psy_lmax = other.psy_lmax;
        this.sfbmax = other.sfbmax;
        this.psymax = other.psymax;
        this.sfbdivide = other.sfbdivide;
        this.width = GrInfo.clone_int(other.width);
        this.window = GrInfo.clone_int(other.window);
        this.count1bits = other.count1bits;

        this.sfb_partition_table = other.sfb_partition_table.slice(0);
        this.slen = GrInfo.clone_int(other.slen);
        this.max_nonzero_coeff = other.max_nonzero_coeff;
    }
}

export { GrInfo };