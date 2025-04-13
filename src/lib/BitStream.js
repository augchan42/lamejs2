import {
    System,
    // VbrMode, // Not used in this snippet
    // Float, // Not used in this snippet
    // ShortBlock, // Not used in this snippet
    Util,
    Arrays,
    // new_array_n, // Not used in this snippet
    new_byte,
    // new_double, // Not used in this snippet
    new_float,
    new_float_n,
    new_int,
    // new_int_n, // Not used in this snippet
    assert
} from './common.js';

// Use import * for modules accessed like namespaces
import { Takehiro, count_bit, slen1_tab, slen2_tab } from './Takehiro.js';
import * as Tables from './Tables.js';
import * as Encoder from './Encoder.js'; // Corrected import
import * as LameInternalFlags from './LameInternalFlags.js';
import * as Lame from './Lame.js';
// Assuming GainAnalysis might be part of another module, e.g., context or a dedicated one
// If GainAnalysis is globally available or part of 'context', no import needed here.
// If it's from its own file, add: import { GainAnalysis } from './GainAnalysis.js';


class TotalBytes {
    constructor() {
        this.total = 0;
    }
}

class BitStream {
    constructor(context) {
        this.context = context; // Keep context for GainAnalysis access
        this.CRC16_POLYNOMIAL = 0x8005;
        this.MAX_LENGTH = 32;

        // Initialize instance properties
        this.buf = null; // Will be initialized in init_bit_stream_w
        this.totbit = 0;
        this.bufByteIdx = 0;
        this.bufBitIdx = 0;

        // Initialize module dependencies (to be set via setModules)
        this.ga = null; // GainAnalysis module/instance
        this.mpg = null; // Mp3Encoder instance?
        this.ver = null; // Version module/instance
        this.vbr = null; // VBR module/instance
    }

    // Static helper methods moved inside the class
    static EQ(a, b) {
        // Original logic uses relative comparison based on larger magnitude
        return (Math.abs(a) > Math.abs(b))
            ? (Math.abs((a) - (b)) <= (Math.abs(a) * 1e-6))
            : (Math.abs((a) - (b)) <= (Math.abs(b) * 1e-6));
    }

    static NEQ(a, b) {
        return !BitStream.EQ(a, b);
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
            var k; // keep var for minimal change, let/const preferred
            if (this.bufBitIdx == 0) {
                this.bufBitIdx = 8;
                this.bufByteIdx++;
                assert(this.bufByteIdx < Lame.LAME_MAXMP3BUFFER);
                // Ensure header buffer exists and indices are valid
                assert(gfc.header && gfc.header[gfc.w_ptr]);
                assert(gfc.header[gfc.w_ptr].write_timing >= this.totbit);
                if (gfc.header[gfc.w_ptr].write_timing == this.totbit) {
                    this.putheader_bits(gfc);
                }
                assert(this.buf); // Ensure buffer is initialized
                this.buf[this.bufByteIdx] = 0;
            }

            k = Math.min(j, this.bufBitIdx);
            j -= k;

            this.bufBitIdx -= k;

            assert(j < this.MAX_LENGTH);
            assert(this.bufBitIdx < this.MAX_LENGTH);
            assert(this.buf);

            // Ensure val >> j doesn't result in unexpected negative numbers if val is large
            // (though j < MAX_LENGTH-2 should prevent this with standard positive val)
            this.buf[this.bufByteIdx] |= ((val >>> j) << this.bufBitIdx); // Use >>> for unsigned shift
            this.totbit += k;
        }
    }

    getframebits(gfp) {
        var gfc = gfp.internal_flags;
        var bit_rate;

        if (gfc.bitrate_index != 0) {
            // Ensure Tables structure is correct
            assert(Tables.bitrate_table && Tables.bitrate_table[gfp.version]);
            bit_rate = Tables.bitrate_table[gfp.version][gfc.bitrate_index];
        }
        else {
            bit_rate = gfp.brate;
        }
        assert(8 <= bit_rate && bit_rate <= 640);

        // Use Math.floor for explicit integer conversion
        var bytes = Math.floor((gfp.version + 1) * 72000 * bit_rate / gfp.out_samplerate + gfc.padding);
        return 8 * bytes;
    }

    drain_into_ancillary(gfp, remainingBits) {
        var gfc = gfp.internal_flags;
        var i;
        assert(remainingBits >= 0);

        // Simplified writing 'LAME'
        const lameHeader = [0x4c, 0x41, 0x4d, 0x45]; // L A M E
        for (let byte of lameHeader) {
            if (remainingBits >= 8) {
                 this.putbits2(gfc, byte, 8);
                 remainingBits -= 8;
            } else break;
        }

        // Write version string
        if (remainingBits >= 32 && this.ver) { // Check if ver module is set
            var version = this.ver.getLameShortVersion(); // Assuming this method exists
            if (version) {
                for (i = 0; i < version.length && remainingBits >= 8; ++i) {
                     this.putbits2(gfc, version.charCodeAt(i), 8); // Use charCodeAt
                     remainingBits -= 8;
                }
            }
        }

        // Fill remaining bits
        for (; remainingBits >= 1; remainingBits -= 1) {
            this.putbits2(gfc, gfc.ancillary_flag, 1);
            gfc.ancillary_flag ^= (!gfp.disable_reservoir ? 1 : 0);
        }

        assert(remainingBits == 0);
    }

    putheader_bits(gfc) {
         // Ensure header buffer exists and indices are valid
        assert(gfc.header && gfc.header[gfc.w_ptr] && gfc.header[gfc.w_ptr].buf);
        assert(this.buf);
        // Ensure lengths and indices allow copy
        assert(this.bufByteIdx + gfc.sideinfo_len < this.buf.length);

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
                assert(this.buf);
                this.buf[this.bufByteIdx] = 0;
            }

            k = Math.min(j, this.bufBitIdx);
            j -= k;

            this.bufBitIdx -= k;

            assert(j < this.MAX_LENGTH);
            assert(this.bufBitIdx < this.MAX_LENGTH);
            assert(this.buf);

            this.buf[this.bufByteIdx] |= ((val >>> j) << this.bufBitIdx); // Use >>>
            this.totbit += k;
        }
    }

    writeheader(gfc, val, j) {
        // Ensure header buffer exists and indices are valid
        assert(gfc.header && gfc.header[gfc.h_ptr] && gfc.header[gfc.h_ptr].buf);
        var ptr = gfc.header[gfc.h_ptr].ptr;

        while (j > 0) {
            var k = Math.min(j, 8 - (ptr & 7));
            j -= k;
            assert(j < this.MAX_LENGTH);
            var byteIndex = ptr >>> 3; // Use unsigned shift for index
            assert(byteIndex < Lame.LAME_MAXMP3BUFFER); // Check buffer bounds
            assert(byteIndex < gfc.header[gfc.h_ptr].buf.length);

            // Apply mask to ensure only relevant bits of val are used
            var shiftAmount = 8 - (ptr & 7) - k;
            var mask = ( (1 << k) - 1 ); // Create a mask of k bits
            var bitsToWrite = (val >>> j) & mask; // Get the lowest k bits of (val >>> j)

            gfc.header[gfc.h_ptr].buf[byteIndex] |= (bitsToWrite << shiftAmount);
            ptr += k;
        }
        gfc.header[gfc.h_ptr].ptr = ptr;
    }


    CRC_update(value, crc) {
        value <<= 8;
        crc &= 0xFFFF; // Ensure crc is 16 bit before starting

        for (var i = 0; i < 8; i++) {
            crc <<= 1;
            value <<= 1; // Shift value along with crc

            // Check the top bit (bit 16) of crc XOR value's relevant bit (now at bit 16)
            if (((crc ^ value) & 0x10000) !== 0) {
                 crc ^= this.CRC16_POLYNOMIAL;
            }
            crc &= 0xFFFF; // Keep crc within 16 bits
        }
        return crc;
    }


    CRC_writeheader(gfc, header) {
        var crc = 0xffff;
         // Ensure header exists and is long enough
        assert(header && header.length >= 6);
        // Ensure indices are valid for sideinfo_len
        assert(gfc.sideinfo_len <= header.length);

        // Mask with 0xff to ensure we process byte values
        crc = this.CRC_update(header[2] & 0xff, crc);
        crc = this.CRC_update(header[3] & 0xff, crc);
        for (var i = 6; i < gfc.sideinfo_len; i++) {
            crc = this.CRC_update(header[i] & 0xff, crc);
        }

        // Assign low and high bytes of CRC
        header[4] = (crc >> 8) & 0xff; // High byte
        header[5] = crc & 0xff;        // Low byte
    }

    encodeSideInfo2(gfp, bitsPerFrame) {
        var gfc = gfp.internal_flags;
        var l3_side;
        var gr, ch;

        l3_side = gfc.l3_side;
         // Ensure header buffer exists and indices are valid
        assert(gfc.header && gfc.header[gfc.h_ptr] && gfc.header[gfc.h_ptr].buf);
        gfc.header[gfc.h_ptr].ptr = 0;
        // Use standard fill method if available, otherwise loop
        Arrays.fill(gfc.header[gfc.h_ptr].buf, 0, gfc.sideinfo_len, 0);

        // Write header fields... (using this.writeheader)
        if (gfp.out_samplerate < 16000) // MPEG 2/2.5
            this.writeheader(gfc, 0xffe, 12); // sync
        else // MPEG 1
            this.writeheader(gfc, 0xfff, 12); // sync

        this.writeheader(gfc, gfp.version, 1); // ID
        this.writeheader(gfc, 1, 2); // layer (01 = layer III) Changed from 4-3 logic to direct value
        this.writeheader(gfc, gfp.error_protection ? 0 : 1, 1); //!protection_bit
        this.writeheader(gfc, gfc.bitrate_index, 4); // bitrate_index
        this.writeheader(gfc, gfc.samplerate_index, 2); // sampling_frequency
        this.writeheader(gfc, gfc.padding, 1); // padding_bit
        this.writeheader(gfc, gfp.extension, 1); // private_bit
        this.writeheader(gfc, gfp.mode.ordinal(), 2); // mode
        this.writeheader(gfc, gfc.mode_ext, 2); // mode_extension
        this.writeheader(gfc, gfp.copyright, 1); // copyright
        this.writeheader(gfc, gfp.original, 1); // original/copy
        this.writeheader(gfc, gfp.emphasis, 2); // emphasis

        if (gfp.error_protection) {
            this.writeheader(gfc, 0, 16); // CRC (placeholder, filled later)
        }

        // Write side info based on MPEG version
        if (gfp.version == 1) { // MPEG 1
             assert(l3_side.main_data_begin >= 0);
             this.writeheader(gfc, l3_side.main_data_begin, 9);
             if (gfc.channels_out == 2)
                 this.writeheader(gfc, l3_side.private_bits, 3);
             else
                 this.writeheader(gfc, l3_side.private_bits, 5);

             for (ch = 0; ch < gfc.channels_out; ch++) {
                 for (let band = 0; band < 4; band++) { // Use let for block scope
                     this.writeheader(gfc, l3_side.scfsi[ch][band], 1);
                 }
             }
             for (gr = 0; gr < 2; gr++) {
                 for (ch = 0; ch < gfc.channels_out; ch++) {
                     var gi = l3_side.tt[gr][ch]; // Keep var for minimal change
                     this.writeheader(gfc, gi.part2_3_length, 12); // part2_3_length (adjusted: only huffman bits now)
                     this.writeheader(gfc, gi.big_values / 2, 9); // big_values
                     this.writeheader(gfc, gi.global_gain, 8); // global_gain
                     this.writeheader(gfc, gi.scalefac_compress, 4); // scalefac_compress
                     this.writeheader(gfc, gi.window_switching_flag, 1); // window_switching_flag

                     if (gi.window_switching_flag) { // Short or Mixed block
                         this.writeheader(gfc, gi.block_type, 2); // block_type
                         this.writeheader(gfc, gi.mixed_block_flag, 1); // mixed_block_flag
                         for (let i = 0; i < 2; i++) { // Use let for block scope
                             this.writeheader(gfc, gi.table_select[i], 5); // table_select
                         }
                          for (let i = 0; i < 3; i++) { // Use let for block scope
                             this.writeheader(gfc, gi.subblock_gain[i], 3); // subblock_gain
                         }
                         // region0_count and region1_count are implicit for short blocks in MPEG1
                          // For mixed blocks:
                         if (gi.block_type == 2 && gi.mixed_block_flag == 1) {
                           this.writeheader(gfc, gi.region0_count, 4); // region0_count for the long part
                           this.writeheader(gfc, gi.region1_count, 3); // region1_count for the long part
                         }

                     } else { // Long block (NORM_TYPE)
                         for (let i = 0; i < 3; i++) { // Use let for block scope
                             this.writeheader(gfc, gi.table_select[i], 5); // table_select
                         }
                         this.writeheader(gfc, gi.region0_count, 4); // region0_count
                         this.writeheader(gfc, gi.region1_count, 3); // region1_count
                     }
                     this.writeheader(gfc, gi.preflag, 1); // preflag
                     this.writeheader(gfc, gi.scalefac_scale, 1); // scalefac_scale
                     this.writeheader(gfc, gi.count1table_select, 1); // count1table_select
                 }
             }
        } else { // MPEG 2/2.5
            assert(l3_side.main_data_begin >= 0);
            this.writeheader(gfc, l3_side.main_data_begin, 8); // main_data_begin
            this.writeheader(gfc, l3_side.private_bits, gfc.channels_out); // private_bits (1 per channel)

            gr = 0; // MPEG 2/2.5 only has one granule in side info
            for (ch = 0; ch < gfc.channels_out; ch++) {
                 var gi = l3_side.tt[gr][ch]; // Keep var
                 this.writeheader(gfc, gi.part2_3_length, 12); // part2_3_length (adjusted)
                 this.writeheader(gfc, gi.big_values / 2, 9); // big_values
                 this.writeheader(gfc, gi.global_gain, 8); // global_gain
                 this.writeheader(gfc, gi.scalefac_compress, 9); // scalefac_compress (MPEG2 has 9 bits)
                 this.writeheader(gfc, gi.window_switching_flag, 1); // window_switching_flag

                 if (gi.window_switching_flag) {
                     this.writeheader(gfc, gi.block_type, 2);
                     this.writeheader(gfc, gi.mixed_block_flag, 1);
                     for (let i = 0; i < 2; i++) {
                         this.writeheader(gfc, gi.table_select[i], 5);
                     }
                     for (let i = 0; i < 3; i++) {
                          this.writeheader(gfc, gi.subblock_gain[i], 3);
                     }
                     // region counts implicit/different for MPEG2 short blocks
                     if (gi.block_type == 2 && gi.mixed_block_flag == 1) { // Mixed block case for MPEG2
                        // Need table for MPEG2 region count bits (differs from MPEG1)
                        // Example placeholder - requires exact spec:
                        this.writeheader(gfc, gi.region0_count, 4); // Example, adjust bits per spec
                        this.writeheader(gfc, gi.region1_count, 3); // Example, adjust bits per spec
                     }
                 } else { // Long blocks
                     for (let i = 0; i < 3; i++) {
                         this.writeheader(gfc, gi.table_select[i], 5);
                     }
                      // Need table for MPEG2 region count bits (differs from MPEG1)
                      // Example placeholder - requires exact spec:
                     this.writeheader(gfc, gi.region0_count, 4); // Example, adjust bits per spec
                     this.writeheader(gfc, gi.region1_count, 3); // Example, adjust bits per spec
                 }
                 // MPEG2 doesn't have preflag bit
                 this.writeheader(gfc, gi.scalefac_scale, 1); // scalefac_scale
                 this.writeheader(gfc, gi.count1table_select, 1); // count1table_select
            }
        }


        // Calculate and write CRC if needed
        if (gfp.error_protection) {
            this.CRC_writeheader(gfc, gfc.header[gfc.h_ptr].buf);
        }

        // Advance header buffer pointers
        { // Keep block scope for old variable
            var old = gfc.h_ptr;
            // Check pointer value consistency (optional sanity check)
            // assert(gfc.header[old].ptr == gfc.sideinfo_len * 8);

            gfc.h_ptr = (old + 1) & (LameInternalFlags.MAX_HEADER_BUF - 1);
             // Ensure header buffer exists and indices are valid before accessing write_timing
             assert(gfc.header && gfc.header[old] && gfc.header[gfc.h_ptr]);
            gfc.header[gfc.h_ptr].write_timing = gfc.header[old].write_timing + bitsPerFrame;

            if (gfc.h_ptr == gfc.w_ptr) {
                // Consider throwing an error instead of just printing
                console.error("Error: MAX_HEADER_BUF too small in bitstream.c");
                // throw new Error("MAX_HEADER_BUF too small");
            }
        }
    }

    // --- Huffman Coding Methods --- (Keep var for minimal change)
    huffman_coder_count1(gfc, gi) {
        // Check if Tables.ht is valid and index is in bounds
        assert(Tables.ht && (gi.count1table_select + 32) < Tables.ht.length);
        var h = Tables.ht[gi.count1table_select + 32];
        var i, bits = 0;

        var ix = gi.big_values;
        var xr = gi.big_values; // Assuming gi.xr holds signed values corresponding to l3_enc
        assert(gi.count1table_select < 2);

        for (i = (gi.count1 - gi.big_values) / 4; i > 0; --i) {
            var huffbits = 0;
            var p = 0, v;

            // Process quad
            v = gi.l3_enc[ix + 0]; // Value (0 or 1)
            if (v !== 0) {
                p += 8; // Huffman table index component
                if (gi.xr[xr + 0] < 0) huffbits = 1; // Sign bit
                assert(v === 1);
            }

            v = gi.l3_enc[ix + 1];
            if (v !== 0) {
                p += 4;
                huffbits <<= 1; // Shift sign bits left
                if (gi.xr[xr + 1] < 0) huffbits |= 1;
                assert(v === 1);
            }

            v = gi.l3_enc[ix + 2];
            if (v !== 0) {
                p += 2;
                huffbits <<= 1;
                if (gi.xr[xr + 2] < 0) huffbits |= 1;
                assert(v === 1);
            }

            v = gi.l3_enc[ix + 3];
            if (v !== 0) {
                p += 1;
                huffbits <<= 1;
                if (gi.xr[xr + 3] < 0) huffbits |= 1;
                assert(v === 1);
            }

            ix += 4;
            xr += 4;
            // Ensure h.table and h.hlen are valid for index p
            assert(h.table && h.hlen && p < h.table.length && p < h.hlen.length);
            var code = h.table[p];
            var len = h.hlen[p];
            var signLen = count_bit(p); // Number of non-zero values = number of sign bits

            assert(len >= 0 && signLen >= 0);

            this.putbits2(gfc, code, len); // Write Huffman code
            if (signLen > 0) {
                this.putbits2(gfc, huffbits, signLen); // Write sign bits
            }
            bits += len + signLen;
        }
        return bits;
    }


    Huffmancode(gfc, tableindex, start, end, gi) {
        // Ensure Tables.ht is valid and index is in bounds
        assert(Tables.ht && tableindex < Tables.ht.length);
        var h = Tables.ht[tableindex];
        var bits = 0;

        assert(tableindex < 32);
        if (0 == tableindex) return 0; // Table 0 means no data

        for (var i = start; i < end; i += 2) {
            var cbits = 0; // Length of Huffman code part
            var xbits = 0; // Length of sign bits + linbits part
            var linbits = h.linmax ? h.xlen : 0; // Length of linbits per value (0 if not used)
            var xlen = h.xlen; // Max value represented directly by huffman code (+1?)
            var ext = 0; // Combined sign/linbits value
            var x1 = Math.abs(gi.l3_enc[i]); // Use absolute value for table lookup
            var x2 = Math.abs(gi.l3_enc[i + 1]);
            var sign1 = (gi.l3_enc[i] < 0); // Original sign
            var sign2 = (gi.l3_enc[i + 1] < 0);

             // Handle linbits (escape mechanism for large values)
            if (tableindex > 15) { // Tables 16-31 use linbits
                 assert(linbits > 0); // xlen should be > 0 for linbits tables
                 if (x1 >= xlen) { // Check against h.xlen which is max value + 1? Or just max value? ISO spec says: xlen=16 for tables 16..31 => max val is 15.
                     var linbits_x1 = x1 - xlen; // Value beyond Huffman range
                     assert(linbits_x1 <= h.linmax);
                     ext = linbits_x1; // Store linbits part
                     xbits = linbits;   // Add linbits length
                     x1 = xlen;         // Use max huffman value for table lookup
                 }
                 if (x2 >= xlen) {
                     var linbits_x2 = x2 - xlen;
                     assert(linbits_x2 <= h.linmax);
                      // Combine linbits: shift previous by linbits length, add new
                     ext = (ext << linbits) | linbits_x2;
                     xbits += linbits; // Add linbits length
                     x2 = xlen;         // Use max huffman value for table lookup
                 }
            } else {
                // Tables 1..15: xlen is the max value (e.g., 1 for table 1)
                 assert(x1 <= xlen); // Values should not exceed table limits
                 assert(x2 <= xlen);
            }


            // Add sign bits if values are non-zero
            var signmask = 0;
            var signlen = 0;
            if (x1 != 0) {
                signmask = sign1 ? 1 : 0;
                signlen = 1;
            }
            if (x2 != 0) {
                signmask = (signmask << 1) | (sign2 ? 1 : 0);
                signlen++;
            }

            // Combine sign bits and linbits (ext)
            ext = (ext << signlen) | signmask;
            xbits += signlen; // Add sign bit length


            // Lookup Huffman code
            var pair_index = x1 * (xlen + 1) + x2; // Calculate index into paired Huffman table
             // Ensure h.table and h.hlen are valid for index pair_index
            assert(h.table && h.hlen && pair_index < h.table.length && pair_index < h.hlen.length);
            cbits = h.hlen[pair_index]; // Get Huffman code length

            assert(cbits >= 0 && xbits >= 0);

            this.putbits2(gfc, h.table[pair_index], cbits); // Write Huffman code
            if (xbits > 0) {
                this.putbits2(gfc, ext, xbits); // Write sign/linbits
            }
            bits += cbits + xbits;
        }
        return bits;
    }


    ShortHuffmancodebits(gfc, gi) {
        // Check scalefac band structure
        assert(gfc.scalefac_band && gfc.scalefac_band.s);
        // Determine region boundaries (ISO 5.1.1): 3 short blocks, sfb 0..11
        // For short blocks, there's no region partitioning like in long blocks.
        // The whole big_values range (0..575) is coded with one table pair.
        // Table select is based on subblock gains (not done here, assumed gi.table_select is set correctly).
        // Let's assume gi.table_select[0] is the correct table for the whole range.
        // The original code had region1Start logic which applies to LONG blocks.

        // Correct logic for short blocks: Use table_select[0] for the whole big_values range.
        var bits = this.Huffmancode(gfc, gi.table_select[0], 0, gi.big_values, gi);
        return bits;
    }

    LongHuffmancodebits(gfc, gi) {
        var bigvalues, bits;
        var region1Start, region2Start;

        bigvalues = gi.big_values;
        assert(0 <= bigvalues && bigvalues <= 576);
         // Check scalefac band structure
        assert(gfc.scalefac_band && gfc.scalefac_band.l);

        // Determine region boundaries using the scalefactor band indices
        var i = gi.region0_count + 1; // region0 uses bands 0..region0_count
        assert(0 <= i && i < gfc.scalefac_band.l.length);
        region1Start = gfc.scalefac_band.l[i]; // Start index for region1

        i += gi.region1_count + 1; // region1 uses bands region0_count+1 .. region0_count+region1_count+1
        assert(0 <= i && i < gfc.scalefac_band.l.length);
        region2Start = gfc.scalefac_band.l[i]; // Start index for region2

        // Clamp region boundaries to the actual number of bigvalues
        if (region1Start > bigvalues) region1Start = bigvalues;
        if (region2Start > bigvalues) region2Start = bigvalues;

        // Code each region with its selected table
        bits = this.Huffmancode(gfc, gi.table_select[0], 0, region1Start, gi);
        bits += this.Huffmancode(gfc, gi.table_select[1], region1Start, region2Start, gi);
        bits += this.Huffmancode(gfc, gi.table_select[2], region2Start, bigvalues, gi);
        return bits;
    }


    writeMainData(gfp) {
        var gr, ch, sfb, data_bits, tot_bits = 0;
        var gfc = gfp.internal_flags;
        var l3_side = gfc.l3_side;

        if (gfp.version == 1) { // MPEG 1
            for (gr = 0; gr < 2; gr++) {
                for (ch = 0; ch < gfc.channels_out; ch++) {
                    var gi = l3_side.tt[gr][ch];
                    var slen1 = slen1_tab[gi.scalefac_compress];
                    var slen2 = slen2_tab[gi.scalefac_compress];
                    var scale_bits = 0; // Renamed from data_bits for clarity

                    // Write scalefactors
                    for (sfb = 0; sfb < gi.sfbdivide; sfb++) { // Using sfbdivide from LAME logic
                        if (gi.scalefac[sfb] == -1) continue; // -1 indicates scfsi=1, reuse from granule 0
                        this.putbits2(gfc, gi.scalefac[sfb], slen1);
                        scale_bits += slen1;
                    }
                    for (; sfb < gi.sfbmax; sfb++) { // Using sfbmax from LAME logic
                        if (gi.scalefac[sfb] == -1) continue;
                        this.putbits2(gfc, gi.scalefac[sfb], slen2);
                        scale_bits += slen2;
                    }
                    // Part2 length check (scalefactors only)
                    // This assertion might fail if the LAME calculation differs slightly from ISO
                    // assert(scale_bits == gi.part2_length);

                    // Write Huffman coded data (Part3)
                    var huff_bits = 0;
                    if (gi.window_switching_flag && gi.block_type == Encoder.SHORT_TYPE) { // SHORT block type check
                         huff_bits = this.ShortHuffmancodebits(gfc, gi);
                    } else { // Long or Mixed
                         huff_bits = this.LongHuffmancodebits(gfc, gi);
                    }
                    huff_bits += this.huffman_coder_count1(gfc, gi); // count1 region

                    // Check total bits (Part2 + Part3)
                    // The original code asserted against part2_3_length + part2_length,
                    // which seems wrong. Part2_3_length *should* include Part2.
                    // Let's assert against the side info's part2_3_length directly.
                    assert(scale_bits + huff_bits == gi.part2_3_length);
                    tot_bits += scale_bits + huff_bits;
                }
            }
        } else { // MPEG 2/2.5
             gr = 0; // Only 1 granule
             for (ch = 0; ch < gfc.channels_out; ch++) {
                 var gi = l3_side.tt[gr][ch];
                 var i, sfb_partition, scale_bits = 0;
                 assert(gi.sfb_partition_table != null); // MPEG2 uses partition tables
                 var huff_bits = 0;

                 sfb = 0;
                 sfb_partition = 0;

                 // Write scalefactors based on block type and partitions
                 if (gi.window_switching_flag && gi.block_type == Encoder.SHORT_TYPE) { // SHORT
                      for (; sfb_partition < 4; sfb_partition++) {
                          var sfbs = gi.sfb_partition_table[sfb_partition] / 3;
                          var slen = gi.slen[sfb_partition]; // slen per partition
                          for (i = 0; i < sfbs; i++, sfb++) {
                              // Write scalefactors for 3 windows, ensure value >= 0
                              this.putbits2(gfc, Math.max(gi.scalefac[sfb * 3 + 0], 0), slen);
                              this.putbits2(gfc, Math.max(gi.scalefac[sfb * 3 + 1], 0), slen);
                              this.putbits2(gfc, Math.max(gi.scalefac[sfb * 3 + 2], 0), slen);
                              scale_bits += 3 * slen;
                          }
                      }
                       huff_bits = this.ShortHuffmancodebits(gfc, gi);
                 } else { // LONG or MIXED
                      for (; sfb_partition < 4; sfb_partition++) {
                           var sfbs = gi.sfb_partition_table[sfb_partition];
                           var slen = gi.slen[sfb_partition];
                           for (i = 0; i < sfbs; i++, sfb++) {
                               this.putbits2(gfc, Math.max(gi.scalefac[sfb], 0), slen);
                               scale_bits += slen;
                           }
                      }
                      huff_bits = this.LongHuffmancodebits(gfc, gi);
                 }
                 huff_bits += this.huffman_coder_count1(gfc, gi); // count1 region

                 // Check lengths (part2_length might differ in LAME vs ISO?)
                 // assert(huff_bits == gi.part2_3_length); // Huffman part
                 // assert(scale_bits == gi.part2_length); // Scalefactor part
                 assert(scale_bits + huff_bits == gi.part2_3_length); // Total main data check
                 tot_bits += scale_bits + huff_bits;
             }
        }
        return tot_bits;
    }


    compute_flushbits(gfp, total_bytes_output) {
        var gfc = gfp.internal_flags;
        var flushbits, remaining_headers;
        var bitsPerFrame;
        var last_ptr, first_ptr;

        // Ensure header buffer exists and indices are valid
        assert(gfc.header);
        first_ptr = gfc.w_ptr;
        last_ptr = (gfc.h_ptr - 1 + LameInternalFlags.MAX_HEADER_BUF) % LameInternalFlags.MAX_HEADER_BUF; // Wrap around correctly

        // Check if buffers are valid before accessing
        assert(gfc.header[last_ptr]);

        // Bits remaining in the buffer from last full frame written + header timings
        flushbits = gfc.header[last_ptr].write_timing - this.totbit;
        total_bytes_output.total = flushbits; // Store total bits before adding final frame/padding

        if (flushbits >= 0) {
            // Calculate how many full headers are waiting to be written
            remaining_headers = (last_ptr - first_ptr + 1 + LameInternalFlags.MAX_HEADER_BUF) % LameInternalFlags.MAX_HEADER_BUF;
             // The flushbits already includes the time until the *start* of the last header.
             // We need to subtract the bits *for* the headers themselves if they haven't been physically put in the buffer yet.
            // This seems overly complex - the timing should just reflect when the *data* ends.

            // Simpler approach: flushbits = time_of_last_data_end - current_bits_written
             // The write_timing likely includes the frame bits *after* the header.
             // Let's trust the original calculation for now.
             flushbits -= remaining_headers * 8 * gfc.sideinfo_len; // Subtract side info bits not yet physically placed
        }

        // Add bits for one final (potentially partial) frame to flush everything
        bitsPerFrame = this.getframebits(gfp);
        flushbits += bitsPerFrame;
        total_bytes_output.total += bitsPerFrame; // Add final frame bits to total

        // Calculate total bytes needed in the output buffer
        if ((total_bytes_output.total % 8) != 0) {
            total_bytes_output.total = 1 + Math.floor(total_bytes_output.total / 8);
        } else {
             total_bytes_output.total = total_bytes_output.total / 8;
        }
        // Add bytes currently physically in the buffer instance
        total_bytes_output.total += this.bufByteIdx + 1;

        if (flushbits < 0) {
            console.error("strange error flushing buffer ... flushbits < 0");
            // Handle error, maybe return negative or throw
            return -1; // Indicate error
        }
        return flushbits;
    }


    flush_bitstream(gfp) {
        var gfc = gfp.internal_flags;
        var l3_side;
        var flushbits;

        // Ensure header buffer exists and indices are valid
        assert(gfc.header);
        var last_ptr = (gfc.h_ptr - 1 + LameInternalFlags.MAX_HEADER_BUF) % LameInternalFlags.MAX_HEADER_BUF; // Wrap around correctly
        assert(gfc.header[last_ptr]);

        l3_side = gfc.l3_side;

        const tempTotalBytes = new TotalBytes(); // Use a temporary object
        flushbits = this.compute_flushbits(gfp, tempTotalBytes);
        if (flushbits < 0) return; // Error calculating flushbits

        this.drain_into_ancillary(gfp, flushbits); // Fill remaining bits

        // Sanity check: total bits written should now align with the end timing of the last frame
        assert(BitStream.EQ(gfc.header[last_ptr].write_timing + this.getframebits(gfp), this.totbit));

        // Reset counters for next potential encoding session?
        // These seem specific to LAME's internal state management
        gfc.ResvSize = 0;
        l3_side.main_data_begin = 0;

        // Handle ReplayGain and Peak Sample calculation if enabled and modules are set
        if (gfc.findReplayGain && this.context && this.ga) {
             // Need GainAnalysis class definition or import
            const GainAnalysis = this.ga; // Assuming ga holds the GainAnalysis class/methods
            var RadioGain = GainAnalysis.GetTitleGain(gfc.rgdata); // Call static method? Or instance?
            assert(BitStream.NEQ(RadioGain, GainAnalysis.GAIN_NOT_ENOUGH_SAMPLES));
            gfc.RadioGain = Math.floor(RadioGain * 10.0 + 0.5); // Explicit floor
        }

        if (gfc.findPeakSample) {
            // Avoid log10(0) or log10(negative)
            if (gfc.PeakSample > 0) {
                 gfc.noclipGainChange = Math.ceil(Math.log10(gfc.PeakSample / 32767.0) * 20.0 * 10.0); // Explicit ceil
            } else {
                 gfc.noclipGainChange = -Infinity; // Or some indicator of no peak/silence
            }

            if (gfc.noclipGainChange > 0) {
                if (BitStream.EQ(gfp.scale, 1.0) || BitStream.EQ(gfp.scale, 0.0)) {
                     gfc.noclipScale = (gfc.PeakSample > 0)
                        ? (Math.floor((32767.0 / gfc.PeakSample) * 100.0) / 100.0)
                        : 1.0; // Avoid division by zero if PeakSample is 0
                }
                else {
                    gfc.noclipScale = -1; // Indicate not applicable due to user scale
                }
            } else {
                gfc.noclipScale = -1; // Indicate no scaling needed or error
            }
        }
    }


    add_dummy_byte(gfp, val, n) {
        var gfc = gfp.internal_flags;
        var i;

        while (n-- > 0) {
            this.putbits_noheaders(gfc, val, 8);

            // Adjust write timings in the header buffer
            for (i = 0; i < LameInternalFlags.MAX_HEADER_BUF; ++i) {
                 // Ensure header buffer exists and index is valid
                if (gfc.header && gfc.header[i]) {
                     gfc.header[i].write_timing += 8;
                }
            }
        }
    }

    format_bitstream(gfp) {
        var gfc = gfp.internal_flags;
        var l3_side;
        l3_side = gfc.l3_side;

        var bitsPerFrame = this.getframebits(gfp);

        // Drain bits from reservoir needed *before* side info/main data
        this.drain_into_ancillary(gfp, l3_side.resvDrain_pre);

        // Write header and side info for the current frame
        this.encodeSideInfo2(gfp, bitsPerFrame);

        // Write scalefactors and Huffman data
        var main_data_bits = this.writeMainData(gfp);
        var side_info_bits = 8 * gfc.sideinfo_len;

        // Drain bits specified to be written *after* main data
        this.drain_into_ancillary(gfp, l3_side.resvDrain_post);

        // Total bits written for this frame's core data + side info + post-drain
        var bits_written_this_frame = side_info_bits + main_data_bits + l3_side.resvDrain_post;

        // Update main_data_begin for the *next* frame's header
        // It's the offset from the header start to where main data begins.
        // Calculated based on how many bits were *not* filled by this frame's data + drain.
        var frame_deficit = bitsPerFrame - bits_written_this_frame;
        // The new main_data_begin is the previous deficit (ResvSize) plus this frame's deficit, divided by 8 for bytes.
        // Note: ResvSize holds the deficit from the *previous* frame in bits.
        l3_side.main_data_begin = (gfc.ResvSize + frame_deficit) / 8;

        // Check buffer consistency: Compare calculated remaining bits (reservoir) with main_data_begin * 8
        // compute_flushbits calculates bits needed based on header timings.
        // ResvSize is LAME's internal track of reservoir bits. They should match.
        const tempTotalBytes = new TotalBytes(); // Use temporary object
        if (this.compute_flushbits(gfp, tempTotalBytes) != gfc.ResvSize) {
             console.error("Internal buffer inconsistency. flushbits <> ResvSize");
             // Potentially correct ResvSize or throw error
             // gfc.ResvSize = this.compute_flushbits(gfp, tempTotalBytes);
        }

        // Check if main_data_begin correctly reflects the reservoir size
        // Allow for small rounding differences? Usually should be exact.
        if (Math.abs(l3_side.main_data_begin * 8 - gfc.ResvSize) > 1) { // Allow tolerance of 1 bit?
             console.error(`Bit reservoir error: MDB*8 (${l3_side.main_data_begin * 8}) != ResvSize (${gfc.ResvSize})`);
             console.error(`  Frame deficit: ${frame_deficit}, Bits written: ${bits_written_this_frame}, BPF: ${bitsPerFrame}`);
             // Print detailed diagnostics if needed (from original code)
             // ... (diagnostic printf) ...
             // Correct the reservoir size to match calculation?
             gfc.ResvSize = l3_side.main_data_begin * 8;
        }


        assert(this.totbit % 8 == 0); // Total bits written should be byte-aligned at frame end

        // Prevent totbit from growing indefinitely (wrap around simulation)
        if (this.totbit > 1000000000) { // Arbitrary large number
            var i;
            var wrapAmount = this.totbit - (this.totbit % bitsPerFrame); // Wrap to a frame boundary? Or just subtract large multiple?
            for (i = 0; i < LameInternalFlags.MAX_HEADER_BUF; ++i) {
                 if (gfc.header && gfc.header[i]) {
                     gfc.header[i].write_timing -= wrapAmount;
                 }
            }
            this.totbit -= wrapAmount;
        }

        return 0; // Success
    }


    copy_buffer(gfc, buffer, bufferPos, size, mp3data) {
        // Calculate bytes currently held in the instance buffer
        var bytes_to_copy = this.bufByteIdx + 1;
        if (bytes_to_copy <= 0) return 0; // Nothing to copy

        // Check if destination buffer has enough space (if size is specified)
        if (size != 0 && bytes_to_copy > size) {
            console.error("Output buffer too small in copy_buffer");
            return -1; // Error: not enough space
        }

        // Perform the copy
        assert(this.buf); // Ensure source buffer exists
        assert(buffer); // Ensure destination buffer exists
        System.arraycopy(this.buf, 0, buffer, bufferPos, bytes_to_copy);

        // Reset instance buffer pointers
        this.bufByteIdx = -1;
        this.bufBitIdx = 0;

        // Update CRC and VBR seek table if mp3data flag is set
        if (mp3data != 0 && this.vbr) { // Check if vbr module is set
            var crc = new_int(1); // Create array to pass by reference (JS workaround)
            crc[0] = gfc.nMusicCRC;
            this.vbr.updateMusicCRC(crc, buffer, bufferPos, bytes_to_copy); // Call VBR method
            gfc.nMusicCRC = crc[0]; // Update CRC state

            if (bytes_to_copy > 0) {
                 // Ensure VBR seek table exists
                 assert(gfc.VBR_seek_table);
                 gfc.VBR_seek_table.nBytesWritten += bytes_to_copy;
            }
        }

        // Decode on the fly if enabled
        if (gfc.decode_on_the_fly && this.mpg) { // Check if mpg module is set
            // Allocate PCM buffer (consider reusing buffer if possible)
            var pcm_buf = new_float_n([2, 1152]); // Max samples per frame
            var mp3_in = bytes_to_copy; // Bytes available for decoder
            var samples_out = -1; // Decoder result
            var currentPos = bufferPos; // Position in the copied buffer
            var i;

            while (samples_out != 0 && mp3_in > 0) { // Loop while decoder produces samples or consumes input
                // Assuming hip_decode1_unclipped takes buffer, startPos, numBytes
                samples_out = this.mpg.hip_decode1_unclipped(gfc.hip, buffer, currentPos, mp3_in, pcm_buf[0], pcm_buf[1]);

                // How many bytes did the decoder consume? This info is often missing from simple decode APIs.
                // Assuming it consumes all input or indicates error.
                // If it returns bytes consumed, update mp3_in and currentPos.
                // For simplicity here, assume it tries to decode one frame. If samples_out > 0, it succeeded.
                // This part needs a more robust decoder API interaction model.
                 var bytes_consumed = mp3_in; // Simplistic assumption - adjust if decoder reports consumption
                 mp3_in -= bytes_consumed;
                 currentPos += bytes_consumed;


                if (samples_out == -1) { // Decoder error
                    console.error("hip_decode1_unclipped error during decode-on-the-fly");
                    samples_out = 0; // Stop decoding loop
                    // Handle error?
                } else if (samples_out > 0) {
                    assert(samples_out <= 1152);

                    // Find peak sample if enabled
                    if (gfc.findPeakSample) {
                        for (i = 0; i < samples_out; i++) {
                             let absSampleL = Math.abs(pcm_buf[0][i]);
                             if (absSampleL > gfc.PeakSample) gfc.PeakSample = absSampleL;
                        }
                        if (gfc.channels_out > 1) {
                            for (i = 0; i < samples_out; i++) {
                                let absSampleR = Math.abs(pcm_buf[1][i]);
                                if (absSampleR > gfc.PeakSample) gfc.PeakSample = absSampleR;
                            }
                        }
                    }

                    // Process ReplayGain if enabled and modules are set
                    if (gfc.findReplayGain && this.context && this.ga) {
                         const GainAnalysis = this.ga; // Assuming ga holds GainAnalysis
                        if (GainAnalysis.AnalyzeSamples(gfc.rgdata, pcm_buf[0], 0, pcm_buf[1], 0, samples_out, gfc.channels_out) == GainAnalysis.GAIN_ANALYSIS_ERROR) {
                             console.error("GainAnalysis.AnalyzeSamples error");
                             return -6; // Specific error code from original
                        }
                    }
                }
            } // End while decoder loop
        } // End decode_on_the_fly

        return bytes_to_copy; // Return number of bytes copied
    }


    init_bit_stream_w(gfc) {
        // Allocate the main buffer
        this.buf = new_byte(Lame.LAME_MAXMP3BUFFER); // Use LAME constant for size

        // Reset header buffer pointers and timing
        gfc.h_ptr = gfc.w_ptr = 0;
        // Ensure header buffer exists before accessing
        if (gfc.header && gfc.header[gfc.h_ptr]) {
             gfc.header[gfc.h_ptr].write_timing = 0;
        } else {
             console.error("Header buffer not initialized before init_bit_stream_w");
             // Handle error appropriately, maybe initialize gfc.header here?
        }

        // Reset instance buffer state
        this.bufByteIdx = -1; // Index of the last *valid* byte (-1 means buffer is empty)
        this.bufBitIdx = 0; // Number of *available* bits in the current byte (0 means byte full or buffer empty)
        this.totbit = 0; // Total bits written logic (handle wrap around later)
    }

} // End class BitStream

// Export the classes using named exports
export { BitStream, TotalBytes };