var common = require('./common.js');
var System = common.System;
var VbrMode = common.VbrMode;
var Float = common.Float;
var ShortBlock = common.ShortBlock;
var Util = common.Util;
var Arrays = common.Arrays;
var new_array_n = common.new_array_n;
var new_byte = common.new_byte;
var new_double = common.new_double;
var new_float = common.new_float;
var new_float_n = common.new_float_n;
var new_int = common.new_int;
var new_int_n = common.new_int_n;
var assert = common.assert;

var Takehiro = require('./Takehiro.js');
var Tables = require('./Tables.js');
var Encoder = require('./Encoder.js');
var LameInternalFlags = require('./LameInternalFlags.js');
var Lame = require('./Lame.js');

class TotalBytes {
    constructor() {
        this.total = 0;
    }
}

class BitStream {
    constructor(context) {
        this.context = context;
        this.CRC16_POLYNOMIAL = 0x8005;
        this.MAX_LENGTH = 32;

        // Initialize instance properties (previously module-level vars)
        this.buf = null;
        this.totbit = 0;
        this.bufByteIdx = 0;
        this.bufBitIdx = 0;
        
        // Initialize module dependencies
        this.ga = null;
        this.mpg = null;
        this.ver = null;
        this.vbr = null;
    }

    setModules(ga, mpg, ver, vbr) {
        this.ga = ga;
        this.mpg = mpg;
        this.ver = ver;
        this.vbr = vbr;
    }

    putbits2(gfc, val, j) {
        assert(j < this.MAX_LENGTH - 2);

        while (j > 0) {
            var k;
            if (this.bufBitIdx == 0) {
                this.bufBitIdx = 8;
                this.bufByteIdx++;
                assert(this.bufByteIdx < Lame.LAME_MAXMP3BUFFER);
                assert(gfc.header[gfc.w_ptr].write_timing >= this.totbit);
                if (gfc.header[gfc.w_ptr].write_timing == this.totbit) {
                    this.putheader_bits(gfc);
                }
                this.buf[this.bufByteIdx] = 0;
            }

            k = Math.min(j, this.bufBitIdx);
            j -= k;

            this.bufBitIdx -= k;

            assert(j < this.MAX_LENGTH);
            assert(this.bufBitIdx < this.MAX_LENGTH);

            this.buf[this.bufByteIdx] |= ((val >> j) << this.bufBitIdx);
            this.totbit += k;
        }
    }

    getframebits(gfp) {
        var gfc = gfp.internal_flags;
        var bit_rate;

        if (gfc.bitrate_index != 0)
            bit_rate = Tables.bitrate_table[gfp.version][gfc.bitrate_index];
        else
            bit_rate = gfp.brate;
        assert(8 <= bit_rate && bit_rate <= 640);

        var bytes = 0 | (gfp.version + 1) * 72000 * bit_rate / gfp.out_samplerate + gfc.padding;
        return 8 * bytes;
    }

    drain_into_ancillary(gfp, remainingBits) {
        var gfc = gfp.internal_flags;
        var i;
        assert(remainingBits >= 0);

        if (remainingBits >= 8) {
            this.putbits2(gfc, 0x4c, 8);
            remainingBits -= 8;
        }
        if (remainingBits >= 8) {
            this.putbits2(gfc, 0x41, 8);
            remainingBits -= 8;
        }
        if (remainingBits >= 8) {
            this.putbits2(gfc, 0x4d, 8);
            remainingBits -= 8;
        }
        if (remainingBits >= 8) {
            this.putbits2(gfc, 0x45, 8);
            remainingBits -= 8;
        }

        if (remainingBits >= 32) {
            var version = this.ver.getLameShortVersion();
            if (remainingBits >= 32)
                for (i = 0; i < version.length && remainingBits >= 8; ++i) {
                    remainingBits -= 8;
                    this.putbits2(gfc, version.charAt(i), 8);
                }
        }

        for (; remainingBits >= 1; remainingBits -= 1) {
            this.putbits2(gfc, gfc.ancillary_flag, 1);
            gfc.ancillary_flag ^= (!gfp.disable_reservoir ? 1 : 0);
        }

        assert(remainingBits == 0);

    }

    putheader_bits(gfc) {
        System.arraycopy(gfc.header[gfc.w_ptr].buf, 0, this.buf, this.bufByteIdx, gfc.sideinfo_len);
        this.bufByteIdx += gfc.sideinfo_len;
        this.totbit += gfc.sideinfo_len * 8;
        gfc.w_ptr = (gfc.w_ptr + 1) & (LameInternalFlags.MAX_HEADER_BUF - 1);
    }

    putbits_noheaders(gfc, val, j) {
        assert(j < this.MAX_LENGTH - 2);

        while (j > 0) {
            var k;
            if (this.bufBitIdx == 0) {
                this.bufBitIdx = 8;
                this.bufByteIdx++;
                assert(this.bufByteIdx < Lame.LAME_MAXMP3BUFFER);
                this.buf[this.bufByteIdx] = 0;
            }

            k = Math.min(j, this.bufBitIdx);
            j -= k;

            this.bufBitIdx -= k;

            assert(j < this.MAX_LENGTH);
            assert(this.bufBitIdx < this.MAX_LENGTH);

            this.buf[this.bufByteIdx] |= ((val >> j) << this.bufBitIdx);
            this.totbit += k;
        }
    }

    writeheader(gfc, val, j) {
        var ptr = gfc.header[gfc.h_ptr].ptr;

        while (j > 0) {
            var k = Math.min(j, 8 - (ptr & 7));
            j -= k;
            assert(j < this.MAX_LENGTH);
            assert((ptr >> 3) < Lame.LAME_MAXMP3BUFFER);

            gfc.header[gfc.h_ptr].buf[ptr >> 3] |= ((val >> j)) << (8 - (ptr & 7) - k);
            ptr += k;
        }
        gfc.header[gfc.h_ptr].ptr = ptr;
    }

    CRC_update(value, crc) {
        value <<= 8;
        for (var i = 0; i < 8; i++) {
            value <<= 1;
            crc <<= 1;

            if ((((crc ^ value) & 0x10000) != 0))
                crc ^= this.CRC16_POLYNOMIAL;
        }
        return crc;
    }

    CRC_writeheader(gfc, header) {
        var crc = 0xffff;
        crc = this.CRC_update(header[2] & 0xff, crc);
        crc = this.CRC_update(header[3] & 0xff, crc);
        for (var i = 6; i < gfc.sideinfo_len; i++) {
            crc = this.CRC_update(header[i] & 0xff, crc);
        }

        header[4] = (byte)(crc >> 8);
        header[5] = (byte)(crc & 255);
    }

    encodeSideInfo2(gfp, bitsPerFrame) {
        var gfc = gfp.internal_flags;
        var l3_side;
        var gr, ch;

        l3_side = gfc.l3_side;
        gfc.header[gfc.h_ptr].ptr = 0;
        Arrays.fill(gfc.header[gfc.h_ptr].buf, 0, gfc.sideinfo_len, 0);
        if (gfp.out_samplerate < 16000)
            this.writeheader(gfc, 0xffe, 12);
        else
            this.writeheader(gfc, 0xfff, 12);
        this.writeheader(gfc, (gfp.version), 1);
        this.writeheader(gfc, 4 - 3, 2);
        this.writeheader(gfc, (!gfp.error_protection ? 1 : 0), 1);
        this.writeheader(gfc, (gfc.bitrate_index), 4);
        this.writeheader(gfc, (gfc.samplerate_index), 2);
        this.writeheader(gfc, (gfc.padding), 1);
        this.writeheader(gfc, (gfp.extension), 1);
        this.writeheader(gfc, (gfp.mode.ordinal()), 2);
        this.writeheader(gfc, (gfc.mode_ext), 2);
        this.writeheader(gfc, (gfp.copyright), 1);
        this.writeheader(gfc, (gfp.original), 1);
        this.writeheader(gfc, (gfp.emphasis), 2);
        if (gfp.error_protection) {
            this.writeheader(gfc, 0, 16);
        }

        if (gfp.version == 1) {
            assert(l3_side.main_data_begin >= 0);
            this.writeheader(gfc, (l3_side.main_data_begin), 9);

            if (gfc.channels_out == 2)
                this.writeheader(gfc, l3_side.private_bits, 3);
            else
                this.writeheader(gfc, l3_side.private_bits, 5);

            for (ch = 0; ch < gfc.channels_out; ch++) {
                var band;
                for (band = 0; band < 4; band++) {
                    this.writeheader(gfc, l3_side.scfsi[ch][band], 1);
                }
            }

            for (gr = 0; gr < 2; gr++) {
                for (ch = 0; ch < gfc.channels_out; ch++) {
                    var gi = l3_side.tt[gr][ch];
                    this.writeheader(gfc, gi.part2_3_length + gi.part2_length, 12);
                    this.writeheader(gfc, gi.big_values / 2, 9);
                    this.writeheader(gfc, gi.global_gain, 8);
                    this.writeheader(gfc, gi.scalefac_compress, 4);

                    if (gi.block_type != Encoder.NORM_TYPE) {
                        this.writeheader(gfc, 1, 1);
                        this.writeheader(gfc, gi.block_type, 2);
                        this.writeheader(gfc, gi.mixed_block_flag, 1);

                        if (gi.table_select[0] == 14)
                            gi.table_select[0] = 16;
                        this.writeheader(gfc, gi.table_select[0], 5);
                        if (gi.table_select[1] == 14)
                            gi.table_select[1] = 16;
                        this.writeheader(gfc, gi.table_select[1], 5);

                        this.writeheader(gfc, gi.subblock_gain[0], 3);
                        this.writeheader(gfc, gi.subblock_gain[1], 3);
                        this.writeheader(gfc, gi.subblock_gain[2], 3);
                    } else {
                        this.writeheader(gfc, 0, 1);
                        this.writeheader(gfc, gi.table_select[0], 5);
                        if (gi.table_select[1] == 14)
                            gi.table_select[1] = 16;
                        this.writeheader(gfc, gi.table_select[1], 5);
                        if (gi.table_select[2] == 14)
                            gi.table_select[2] = 16;
                        this.writeheader(gfc, gi.table_select[2], 5);

                        assert(0 <= gi.region0_count && gi.region0_count < 16);
                        assert(0 <= gi.region1_count && gi.region1_count < 8);
                        this.writeheader(gfc, gi.region0_count, 4);
                        this.writeheader(gfc, gi.region1_count, 3);
                    }
                    this.writeheader(gfc, gi.preflag, 1);
                    this.writeheader(gfc, gi.scalefac_scale, 1);
                    this.writeheader(gfc, gi.count1table_select, 1);
                }
            }
        } else {
            assert(l3_side.main_data_begin >= 0);
            this.writeheader(gfc, (l3_side.main_data_begin), 8);
            this.writeheader(gfc, l3_side.private_bits, gfc.channels_out);

            gr = 0;
            for (ch = 0; ch < gfc.channels_out; ch++) {
                var gi = l3_side.tt[gr][ch];
                this.writeheader(gfc, gi.part2_3_length + gi.part2_length, 12);
                this.writeheader(gfc, gi.big_values / 2, 9);
                this.writeheader(gfc, gi.global_gain, 8);
                this.writeheader(gfc, gi.scalefac_compress, 9);

                if (gi.block_type != Encoder.NORM_TYPE) {
                    this.writeheader(gfc, 1, 1);
                    this.writeheader(gfc, gi.block_type, 2);
                    this.writeheader(gfc, gi.mixed_block_flag, 1);

                    if (gi.table_select[0] == 14)
                        gi.table_select[0] = 16;
                    this.writeheader(gfc, gi.table_select[0], 5);
                    if (gi.table_select[1] == 14)
                        gi.table_select[1] = 16;
                    this.writeheader(gfc, gi.table_select[1], 5);

                    this.writeheader(gfc, gi.subblock_gain[0], 3);
                    this.writeheader(gfc, gi.subblock_gain[1], 3);
                    this.writeheader(gfc, gi.subblock_gain[2], 3);
                } else {
                    this.writeheader(gfc, 0, 1);
                    this.writeheader(gfc, gi.table_select[0], 5);
                    if (gi.table_select[1] == 14)
                        gi.table_select[1] = 16;
                    this.writeheader(gfc, gi.table_select[1], 5);
                    if (gi.table_select[2] == 14)
                        gi.table_select[2] = 16;
                    this.writeheader(gfc, gi.table_select[2], 5);

                    assert(0 <= gi.region0_count && gi.region0_count < 16);
                    assert(0 <= gi.region1_count && gi.region1_count < 8);
                    this.writeheader(gfc, gi.region0_count, 4);
                    this.writeheader(gfc, gi.region1_count, 3);
                }

                this.writeheader(gfc, gi.scalefac_scale, 1);
                this.writeheader(gfc, gi.count1table_select, 1);
            }
        }

        if (gfp.error_protection) {
            this.CRC_writeheader(gfc, gfc.header[gfc.h_ptr].buf);
        }

        {
            var old = gfc.h_ptr;
            assert(gfc.header[old].ptr == gfc.sideinfo_len * 8);

            gfc.h_ptr = (old + 1) & (LameInternalFlags.MAX_HEADER_BUF - 1);
            gfc.header[gfc.h_ptr].write_timing = gfc.header[old].write_timing
                + bitsPerFrame;

            if (gfc.h_ptr == gfc.w_ptr) {
                System.err
                    .println("Error: MAX_HEADER_BUF too small in bitstream.c \n");
            }

        }
    }

    huffman_coder_count1(gfc, gi) {
        var h = Tables.ht[gi.count1table_select + 32];
        var i, bits = 0;

        var ix = gi.big_values;
        var xr = gi.big_values;
        assert(gi.count1table_select < 2);

        for (i = (gi.count1 - gi.big_values) / 4; i > 0; --i) {
            var huffbits = 0;
            var p = 0, v;

            v = gi.l3_enc[ix + 0];
            if (v != 0) {
                p += 8;
                if (gi.xr[xr + 0] < 0)
                    huffbits++;
                assert(v <= 1);
            }

            v = gi.l3_enc[ix + 1];
            if (v != 0) {
                p += 4;
                huffbits *= 2;
                if (gi.xr[xr + 1] < 0)
                    huffbits++;
                assert(v <= 1);
            }

            v = gi.l3_enc[ix + 2];
            if (v != 0) {
                p += 2;
                huffbits *= 2;
                if (gi.xr[xr + 2] < 0)
                    huffbits++;
                assert(v <= 1);
            }

            v = gi.l3_enc[ix + 3];
            if (v != 0) {
                p++;
                huffbits *= 2;
                if (gi.xr[xr + 3] < 0)
                    huffbits++;
                assert(v <= 1);
            }

            ix += 4;
            xr += 4;
            this.putbits2(gfc, huffbits + h.table[p], h.hlen[p]);
            bits += h.hlen[p];
        }
        return bits;
    }

    Huffmancode(gfc, tableindex, start, end, gi) {
        var h = Tables.ht[tableindex];
        var bits = 0;

        assert(tableindex < 32);
        if (0 == tableindex)
            return bits;

        for (var i = start; i < end; i += 2) {
            var cbits = 0;
            var xbits = 0;
            var linbits = h.xlen;
            var xlen = h.xlen;
            var ext = 0;
            var x1 = gi.l3_enc[i];
            var x2 = gi.l3_enc[i + 1];

            if (x1 != 0) {
                if (gi.xr[i] < 0)
                    ext++;
                cbits--;
            }

            if (tableindex > 15) {
                if (x1 > 14) {
                    var linbits_x1 = x1 - 15;
                    assert(linbits_x1 <= h.linmax);
                    ext |= linbits_x1 << 1;
                    xbits = linbits;
                    x1 = 15;
                }

                if (x2 > 14) {
                    var linbits_x2 = x2 - 15;
                    assert(linbits_x2 <= h.linmax);
                    ext <<= linbits;
                    ext |= linbits_x2;
                    xbits += linbits;
                    x2 = 15;
                }
                xlen = 16;
            }

            if (x2 != 0) {
                ext <<= 1;
                if (gi.xr[i + 1] < 0)
                    ext++;
                cbits--;
            }

            assert((x1 | x2) < 16);

            x1 = x1 * xlen + x2;
            xbits -= cbits;
            cbits += h.hlen[x1];

            assert(cbits <= this.MAX_LENGTH);
            assert(xbits <= this.MAX_LENGTH);

            this.putbits2(gfc, h.table[x1], cbits);
            this.putbits2(gfc, ext, xbits);
            bits += cbits + xbits;
        }
        return bits;
    }

    ShortHuffmancodebits(gfc, gi) {
        var region1Start = 3 * gfc.scalefac_band.s[3];
        if (region1Start > gi.big_values)
            region1Start = gi.big_values;

        var bits = this.Huffmancode(gfc, gi.table_select[0], 0, region1Start, gi);
        bits += this.Huffmancode(gfc, gi.table_select[1], region1Start,
            gi.big_values, gi);
        return bits;
    }

    LongHuffmancodebits(gfc, gi) {
        var bigvalues, bits;
        var region1Start, region2Start;

        bigvalues = gi.big_values;
        assert(0 <= bigvalues && bigvalues <= 576);

        var i = gi.region0_count + 1;
        assert(0 <= i);
        assert(i < gfc.scalefac_band.l.length);
        region1Start = gfc.scalefac_band.l[i];
        i += gi.region1_count + 1;
        assert(0 <= i);
        assert(i < gfc.scalefac_band.l.length);
        region2Start = gfc.scalefac_band.l[i];

        if (region1Start > bigvalues)
            region1Start = bigvalues;

        if (region2Start > bigvalues)
            region2Start = bigvalues;

        bits = this.Huffmancode(gfc, gi.table_select[0], 0, region1Start, gi);
        bits += this.Huffmancode(gfc, gi.table_select[1], region1Start,
            region2Start, gi);
        bits += this.Huffmancode(gfc, gi.table_select[2], region2Start, bigvalues,
            gi);
        return bits;
    }

    writeMainData(gfp) {
        var gr, ch, sfb, data_bits, tot_bits = 0;
        var gfc = gfp.internal_flags;
        var l3_side = gfc.l3_side;

        if (gfp.version == 1) {
            for (gr = 0; gr < 2; gr++) {
                for (ch = 0; ch < gfc.channels_out; ch++) {
                    var gi = l3_side.tt[gr][ch];
                    var slen1 = Takehiro.slen1_tab[gi.scalefac_compress];
                    var slen2 = Takehiro.slen2_tab[gi.scalefac_compress];
                    data_bits = 0;
                    for (sfb = 0; sfb < gi.sfbdivide; sfb++) {
                        if (gi.scalefac[sfb] == -1)
                            continue;
                        this.putbits2(gfc, gi.scalefac[sfb], slen1);
                        data_bits += slen1;
                    }
                    for (; sfb < gi.sfbmax; sfb++) {
                        if (gi.scalefac[sfb] == -1)
                            continue;
                        this.putbits2(gfc, gi.scalefac[sfb], slen2);
                        data_bits += slen2;
                    }
                    assert(data_bits == gi.part2_length);

                    if (gi.block_type == Encoder.SHORT_TYPE) {
                        data_bits += this.ShortHuffmancodebits(gfc, gi);
                    } else {
                        data_bits += this.LongHuffmancodebits(gfc, gi);
                    }
                    data_bits += this.huffman_coder_count1(gfc, gi);
                    assert(data_bits == gi.part2_3_length + gi.part2_length);
                    tot_bits += data_bits;
                }
            }
        } else {
            gr = 0;
            for (ch = 0; ch < gfc.channels_out; ch++) {
                var gi = l3_side.tt[gr][ch];
                var i, sfb_partition, scale_bits = 0;
                assert(gi.sfb_partition_table != null);
                data_bits = 0;
                sfb = 0;
                sfb_partition = 0;

                if (gi.block_type == Encoder.SHORT_TYPE) {
                    for (; sfb_partition < 4; sfb_partition++) {
                        var sfbs = gi.sfb_partition_table[sfb_partition] / 3;
                        var slen = gi.slen[sfb_partition];
                        for (i = 0; i < sfbs; i++, sfb++) {
                            this.putbits2(gfc,
                                Math.max(gi.scalefac[sfb * 3 + 0], 0), slen);
                            this.putbits2(gfc,
                                Math.max(gi.scalefac[sfb * 3 + 1], 0), slen);
                            this.putbits2(gfc,
                                Math.max(gi.scalefac[sfb * 3 + 2], 0), slen);
                            scale_bits += 3 * slen;
                        }
                    }
                    data_bits += this.ShortHuffmancodebits(gfc, gi);
                } else {
                    for (; sfb_partition < 4; sfb_partition++) {
                        var sfbs = gi.sfb_partition_table[sfb_partition];
                        var slen = gi.slen[sfb_partition];
                        for (i = 0; i < sfbs; i++, sfb++) {
                            this.putbits2(gfc, Math.max(gi.scalefac[sfb], 0), slen);
                            scale_bits += slen;
                        }
                    }
                    data_bits += this.LongHuffmancodebits(gfc, gi);
                }
                data_bits += this.huffman_coder_count1(gfc, gi);
                assert(data_bits == gi.part2_3_length);
                assert(scale_bits == gi.part2_length);
                tot_bits += scale_bits + data_bits;
            }
        }
        return tot_bits;
    }

    compute_flushbits(gfp, total_bytes_output) {
        var gfc = gfp.internal_flags;
        var flushbits, remaining_headers;
        var bitsPerFrame;
        var last_ptr, first_ptr;
        first_ptr = gfc.w_ptr;
        last_ptr = gfc.h_ptr - 1;
        if (last_ptr == -1)
            last_ptr = LameInternalFlags.MAX_HEADER_BUF - 1;

        flushbits = gfc.header[last_ptr].write_timing - this.totbit;
        total_bytes_output.total = flushbits;

        if (flushbits >= 0) {
            remaining_headers = 1 + last_ptr - first_ptr;
            if (last_ptr < first_ptr)
                remaining_headers = 1 + last_ptr - first_ptr
                    + LameInternalFlags.MAX_HEADER_BUF;
            flushbits -= remaining_headers * 8 * gfc.sideinfo_len;
        }

        bitsPerFrame = this.getframebits(gfp);
        flushbits += bitsPerFrame;
        total_bytes_output.total += bitsPerFrame;
        if ((total_bytes_output.total % 8) != 0)
            total_bytes_output.total = 1 + (total_bytes_output.total / 8);
        else
            total_bytes_output.total = (total_bytes_output.total / 8);
        total_bytes_output.total += this.bufByteIdx + 1;

        if (flushbits < 0) {
            System.err.println("strange error flushing buffer ... \n");
        }
        return flushbits;
    }

    flush_bitstream(gfp) {
        var gfc = gfp.internal_flags;
        var l3_side;
        var flushbits;
        var last_ptr = gfc.h_ptr - 1;
        if (last_ptr == -1)
            last_ptr = LameInternalFlags.MAX_HEADER_BUF - 1;
        l3_side = gfc.l3_side;

        if ((flushbits = this.compute_flushbits(gfp, new TotalBytes())) < 0)
            return;
        this.drain_into_ancillary(gfp, flushbits);

        assert(gfc.header[last_ptr].write_timing + this.getframebits(gfp) == this.totbit);

        gfc.ResvSize = 0;
        l3_side.main_data_begin = 0;

        if (gfc.findReplayGain) {
            var RadioGain = this.context.gainAnalysis.GetTitleGain(gfc.rgdata);
            assert(this.NEQ(RadioGain, GainAnalysis.GAIN_NOT_ENOUGH_SAMPLES));
            gfc.RadioGain = Math.floor(RadioGain * 10.0 + 0.5) | 0;
        }

        if (gfc.findPeakSample) {
            gfc.noclipGainChange = Math.ceil(Math
                        .log10(gfc.PeakSample / 32767.0) * 20.0 * 10.0) | 0;

            if (gfc.noclipGainChange > 0) {
                if (this.EQ(gfp.scale, 1.0) || this.EQ(gfp.scale, 0.0))
                    gfc.noclipScale = (Math
                        .floor((32767.0 / gfc.PeakSample) * 100.0) / 100.0);
                else {
                    gfc.noclipScale = -1;
                }
            } else
                gfc.noclipScale = -1;
        }
    }

    add_dummy_byte(gfp, val, n) {
        var gfc = gfp.internal_flags;
        var i;

        while (n-- > 0) {
            this.putbits_noheaders(gfc, val, 8);

            for (i = 0; i < LameInternalFlags.MAX_HEADER_BUF; ++i)
                gfc.header[i].write_timing += 8;
        }
    }

    format_bitstream(gfp) {
        var gfc = gfp.internal_flags;
        var l3_side;
        l3_side = gfc.l3_side;

        var bitsPerFrame = this.getframebits(gfp);
        this.drain_into_ancillary(gfp, l3_side.resvDrain_pre);

        this.encodeSideInfo2(gfp, bitsPerFrame);
        var bits = 8 * gfc.sideinfo_len;
        bits += this.writeMainData(gfp);
        this.drain_into_ancillary(gfp, l3_side.resvDrain_post);
        bits += l3_side.resvDrain_post;

        l3_side.main_data_begin += (bitsPerFrame - bits) / 8;

        if (this.compute_flushbits(gfp, new TotalBytes()) != gfc.ResvSize) {
            System.err.println("Internal buffer inconsistency. flushbits <> ResvSize");
        }

        if ((l3_side.main_data_begin * 8) != gfc.ResvSize) {
            System.err.printf("bit reservoir error: \n"
                + "l3_side.main_data_begin: %d \n"
                + "Resvoir size:             %d \n"
                + "resv drain (post)         %d \n"
                + "resv drain (pre)          %d \n"
                + "header and sideinfo:      %d \n"
                + "data bits:                %d \n"
                + "total bits:               %d (remainder: %d) \n"
                + "bitsperframe:             %d \n",
                8 * l3_side.main_data_begin, gfc.ResvSize,
                l3_side.resvDrain_post, l3_side.resvDrain_pre,
                8 * gfc.sideinfo_len, bits - l3_side.resvDrain_post - 8
                * gfc.sideinfo_len, bits, bits % 8, bitsPerFrame);

            System.err.println("This is a fatal error.  It has several possible causes:");
            System.err.println("90%%  LAME compiled with buggy version of gcc using advanced optimizations");
            System.err.println(" 9%%  Your system is overclocked");
            System.err.println(" 1%%  bug in LAME encoding library");

            gfc.ResvSize = l3_side.main_data_begin * 8;
        }
        assert(this.totbit % 8 == 0);

        if (this.totbit > 1000000000) {
            var i;
            for (i = 0; i < LameInternalFlags.MAX_HEADER_BUF; ++i)
                gfc.header[i].write_timing -= this.totbit;
            this.totbit = 0;
        }

        return 0;
    }

    copy_buffer(gfc, buffer, bufferPos, size, mp3data) {
        var minimum = this.bufByteIdx + 1;
        if (minimum <= 0)
            return 0;
        if (size != 0 && minimum > size) {
            return -1;
        }
        System.arraycopy(this.buf, 0, buffer, bufferPos, minimum);
        this.bufByteIdx = -1;
        this.bufBitIdx = 0;

        if (mp3data != 0) {
            var crc = new_int(1);
            crc[0] = gfc.nMusicCRC;
            this.vbr.updateMusicCRC(crc, buffer, bufferPos, minimum);
            gfc.nMusicCRC = crc[0];

            if (minimum > 0) {
                gfc.VBR_seek_table.nBytesWritten += minimum;
            }

            if (gfc.decode_on_the_fly) {
                var pcm_buf = new_float_n([2, 1152]);
                var mp3_in = minimum;
                var samples_out = -1;
                var i;

                while (samples_out != 0) {

                    samples_out = mpg.hip_decode1_unclipped(gfc.hip, buffer,
                        bufferPos, mp3_in, pcm_buf[0], pcm_buf[1]);

                    mp3_in = 0;

                    if (samples_out == -1) {
                        samples_out = 0;
                    }
                    if (samples_out > 0) {
                        assert(samples_out <= 1152);

                        if (gfc.findPeakSample) {
                            for (i = 0; i < samples_out; i++) {
                                if (pcm_buf[0][i] > gfc.PeakSample)
                                    gfc.PeakSample = pcm_buf[0][i];
                                else if (-pcm_buf[0][i] > gfc.PeakSample)
                                    gfc.PeakSample = -pcm_buf[0][i];
                            }
                            if (gfc.channels_out > 1)
                                for (i = 0; i < samples_out; i++) {
                                    if (pcm_buf[1][i] > gfc.PeakSample)
                                        gfc.PeakSample = pcm_buf[1][i];
                                    else if (-pcm_buf[1][i] > gfc.PeakSample)
                                        gfc.PeakSample = -pcm_buf[1][i];
                                }
                        }

                        if (gfc.findReplayGain)
                            if (ga.AnalyzeSamples(gfc.rgdata, pcm_buf[0], 0,
                                    pcm_buf[1], 0, samples_out,
                                    gfc.channels_out) == GainAnalysis.GAIN_ANALYSIS_ERROR)
                                return -6;

                    }
                }
            }
        }
        return minimum;
    }

    init_bit_stream_w(gfc) {
        this.buf = new_byte(Lame.LAME_MAXMP3BUFFER);

        gfc.h_ptr = gfc.w_ptr = 0;
        gfc.header[gfc.h_ptr].write_timing = 0;
        this.bufByteIdx = -1;
        this.bufBitIdx = 0;
        this.totbit = 0;
    }
}

BitStream.EQ = function (a, b) {
    return (Math.abs(a) > Math.abs(b)) 
        ? (Math.abs((a) - (b)) <= (Math.abs(a) * 1e-6))
        : (Math.abs((a) - (b)) <= (Math.abs(b) * 1e-6));
};

BitStream.NEQ = function (a, b) {
    return !BitStream.EQ(a, b);
};

module.exports = BitStream;
module.exports.TotalBytes = TotalBytes;
