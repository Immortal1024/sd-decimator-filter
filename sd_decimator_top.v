// Synthesizable single-clock CIC + FIR decimator.
//
// Architecture:
//   1-bit input -> (-1)^n demodulation -> CIC (N=4, R=512, M=1)
//   -> 17-tap symmetric FIR -> decimate by 4 -> rounded/saturated signed 20-bit output.
//
// Reference rates assumed by the MATLAB model: 4.096 MHz -> 8 kHz -> 2 kHz.
// The module itself is clock-rate agnostic; one input sample is accepted per clk
// when in_valid=1. Use a 4.096 MHz sample clock if the source vectors represent
// 4.096 Msamples/s.
//
// CIC arithmetic wraps in 38-bit two's-complement registers. FIR coefficients
// are the 17-tap weighted-least-squares values from the supplied MATLAB design,
// quantized to signed Q24 integers. Output scaling approximates
// G = 0.90*2^19/0.285 = 1655646.315789... using integer 1655646.

module sd_decimator_top (
    input  wire                  clk,
    input  wire                  reset,       // synchronous active-high reset
    input  wire                  in_valid,
    input  wire                  bit_in,      // unsigned 1-bit input sample
    output reg                   out_valid,   // one-cycle pulse per output sample
    output reg signed [19:0]     out_data
);

    // ---------------------------------------------------------------------
    // Input demodulation: x[n] = (-1)^n * bit_in, with n=0 even.
    // The 1-bit input itself is unsigned; the demodulated value is signed 2-bit.
    // ---------------------------------------------------------------------
    reg input_parity;
    wire signed [1:0] x_demod;
    wire signed [37:0] x_ext;

    assign x_demod = !bit_in ? 2'sd0 : (input_parity ? -2'sd1 : 2'sd1);
    assign x_ext = {{36{x_demod[1]}}, x_demod};

    // ---------------------------------------------------------------------
    // CIC integrators: four 38-bit stages at the input sample rate.
    // Intermediate values intentionally wrap modulo 2^38.
    // ---------------------------------------------------------------------
    reg signed [37:0] integ1, integ2, integ3, integ4;
    wire signed [37:0] integ1_next = integ1 + x_ext;
    wire signed [37:0] integ2_next = integ2 + integ1_next;
    wire signed [37:0] integ3_next = integ3 + integ2_next;
    wire signed [37:0] integ4_next = integ4 + integ3_next;

    // Every 512 accepted input samples, take the newly updated final integrator.
    reg [8:0] cic_count;
    wire cic_tick = in_valid && (cic_count == 9'd511);

    // ---------------------------------------------------------------------
    // CIC combs: four stages, differential delay M=1, run only at CIC rate.
    // Use *_now signals so the current sample propagates through all four combs
    // in this same accepted decimation event (matching the MATLAB filter output).
    // ---------------------------------------------------------------------
    reg signed [37:0] comb_dly0, comb_dly1, comb_dly2, comb_dly3;
    wire signed [37:0] comb0_now = integ4_next - comb_dly0;
    wire signed [37:0] comb1_now = comb0_now - comb_dly1;
    wire signed [37:0] comb2_now = comb1_now - comb_dly2;
    wire signed [37:0] comb3_now = comb2_now - comb_dly3;
    wire signed [37:0] cic_sample_now = comb3_now;

    // ---------------------------------------------------------------------
    // FIR coefficients. Q24 integers rounded from MATLAB's hfir vector.
    // The coefficients are symmetric, allowing 17 taps to be implemented with
    // 8 pre-adders + 8 multipliers + 1 centre multiplier.
    // ---------------------------------------------------------------------
    localparam signed [23:0] C0 = -24'sd59144;
    localparam signed [23:0] C1 = -24'sd29620;
    localparam signed [23:0] C2 =  24'sd296960;
    localparam signed [23:0] C3 =  24'sd347658;
    localparam signed [23:0] C4 = -24'sd676079;
    localparam signed [23:0] C5 = -24'sd1511568;
    localparam signed [23:0] C6 =  24'sd629633;
    localparam signed [23:0] C7 =  24'sd5392870;
    localparam signed [23:0] C8 =  24'sd8009677;

    // History contains the previous 16 CIC samples. hist[0] is x[n-1].
    reg signed [37:0] fir_hist [0:15];

    // Pair symmetric samples. Pair sums require 39 signed bits.
    wire signed [38:0] pair0 = {cic_sample_now[37], cic_sample_now} + {fir_hist[15][37], fir_hist[15]};
    wire signed [38:0] pair1 = {fir_hist[0][37],  fir_hist[0]}  + {fir_hist[14][37], fir_hist[14]};
    wire signed [38:0] pair2 = {fir_hist[1][37],  fir_hist[1]}  + {fir_hist[13][37], fir_hist[13]};
    wire signed [38:0] pair3 = {fir_hist[2][37],  fir_hist[2]}  + {fir_hist[12][37], fir_hist[12]};
    wire signed [38:0] pair4 = {fir_hist[3][37],  fir_hist[3]}  + {fir_hist[11][37], fir_hist[11]};
    wire signed [38:0] pair5 = {fir_hist[4][37],  fir_hist[4]}  + {fir_hist[10][37], fir_hist[10]};
    wire signed [38:0] pair6 = {fir_hist[5][37],  fir_hist[5]}  + {fir_hist[9][37],  fir_hist[9]};
    wire signed [38:0] pair7 = {fir_hist[6][37],  fir_hist[6]}  + {fir_hist[8][37],  fir_hist[8]};

    // Full-width signed products: 39x24 -> 63 bits; centre is 38x24 -> 62 bits.
    wire signed [62:0] p0 = $signed(pair0) * $signed(C0);
    wire signed [62:0] p1 = $signed(pair1) * $signed(C1);
    wire signed [62:0] p2 = $signed(pair2) * $signed(C2);
    wire signed [62:0] p3 = $signed(pair3) * $signed(C3);
    wire signed [62:0] p4 = $signed(pair4) * $signed(C4);
    wire signed [62:0] p5 = $signed(pair5) * $signed(C5);
    wire signed [62:0] p6 = $signed(pair6) * $signed(C6);
    wire signed [62:0] p7 = $signed(pair7) * $signed(C7);
    wire signed [61:0] p8 = $signed(fir_hist[7]) * $signed(C8);

    // 9 terms require 68 bits after sign extension.
    wire signed [67:0] fir_acc_comb =
          {{5{p0[62]}}, p0}
        + {{5{p1[62]}}, p1}
        + {{5{p2[62]}}, p2}
        + {{5{p3[62]}}, p3}
        + {{5{p4[62]}}, p4}
        + {{5{p5[62]}}, p5}
        + {{5{p6[62]}}, p6}
        + {{5{p7[62]}}, p7}
        + {{6{p8[61]}}, p8};

    // Restore the MATLAB normalization and scale to the signed 20-bit range.
    // fir_acc_comb represents filtered CIC raw counts times 2^24.
    // cic raw counts are 2^36 larger than cicn, so divide by 2^(24+36)=2^60.
    localparam signed [21:0] G_INT = 22'sd1655646;
    wire signed [89:0] scaled_product = fir_acc_comb * G_INT;
    localparam signed [89:0] ROUND_HALF = (90'sd1 << 59);
    wire signed [89:0] rounded_value = (scaled_product >= 0)
        ? ((scaled_product + ROUND_HALF) >>> 60)
        : -(((-scaled_product + ROUND_HALF) >>> 60));

    localparam signed [89:0] MAX20_EXT =  90'sd524287;
    localparam signed [89:0] MIN20_EXT = -90'sd524288;
    wire signed [19:0] out_data_next =
        (rounded_value > MAX20_EXT) ? 20'sh7ffff :
        (rounded_value < MIN20_EXT) ? 20'sh80000 :
        rounded_value[19:0];

    // FIR output decimation phase: MATLAB uses yf(4:4:end), i.e. the 4th,
    // 8th, 12th, ... FIR output sample (zero-based phases 3, 7, 11, ...).
    reg [1:0] fir_decim_phase;
    integer j;

    always @(posedge clk) begin
        if (reset) begin
            input_parity <= 1'b0;
            integ1 <= 38'sd0;
            integ2 <= 38'sd0;
            integ3 <= 38'sd0;
            integ4 <= 38'sd0;
            cic_count <= 9'd0;
            comb_dly0 <= 38'sd0;
            comb_dly1 <= 38'sd0;
            comb_dly2 <= 38'sd0;
            comb_dly3 <= 38'sd0;
            fir_decim_phase <= 2'd0;
            out_valid <= 1'b0;
            out_data <= 20'sd0;
            for (j = 0; j < 16; j = j + 1)
                fir_hist[j] <= 38'sd0;
        end else begin
            out_valid <= 1'b0;

            if (in_valid) begin
                // Input parity advances only for accepted samples.
                input_parity <= ~input_parity;

                // Update the four integrators with this input sample.
                integ1 <= integ1_next;
                integ2 <= integ2_next;
                integ3 <= integ3_next;
                integ4 <= integ4_next;

                if (cic_count == 9'd511) begin
                    cic_count <= 9'd0;

                    // Comb delay lines update once per decimated CIC sample.
                    comb_dly0 <= integ4_next;
                    comb_dly1 <= comb0_now;
                    comb_dly2 <= comb1_now;
                    comb_dly3 <= comb2_now;

                    // Feed this current CIC sample into the FIR and shift history.
                    for (j = 15; j > 0; j = j - 1)
                        fir_hist[j] <= fir_hist[j-1];
                    fir_hist[0] <= cic_sample_now;

                    if (fir_decim_phase == 2'd3) begin
                        out_data <= out_data_next;
                        out_valid <= 1'b1;
                        fir_decim_phase <= 2'd0;
                    end else begin
                        fir_decim_phase <= fir_decim_phase + 2'd1;
                    end
                end else begin
                    cic_count <= cic_count + 9'd1;
                end
            end
        end
    end

endmodule
