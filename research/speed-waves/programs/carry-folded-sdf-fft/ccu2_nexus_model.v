// ccu2_nexus_model.v (as-5): the Nexus CCU2 carry cell with the semantics of Yosys lattice/cells_sim_nexus.v
// (S = LUT4 ^ cin; cout = LUT4 ? cin : LUT3(A,B,C) = INIT[7:0]), X-optimistic mux trees, INIT given as "0xHHHH".
module lut4n(input [15:0] I, input A, B, C, D, output Z, output Z3);
  wire [7:0] s3 = D ? I[15:8] : I[7:0]; wire [3:0] s2 = C ? s3[7:4] : s3[3:0];
  wire [1:0] s1 = B ? s2[3:2] : s2[1:0]; assign Z = A ? s1[1] : s1[0];
  wire [3:0] t2 = C ? I[7:4] : I[3:0]; wire [1:0] t1 = B ? t2[3:2] : t2[1:0]; assign Z3 = A ? t1[1] : t1[0];
endmodule
module CCU2(input CIN, input A0, B0, C0, D0, A1, B1, C1, D1, output S0, S1, output COUT);
  parameter INIT0 = "0x0000"; parameter INIT1 = "0x0000"; parameter INJECT = "YES";
  function [15:0] hx(input [47:0] s);
    integer i; reg [7:0] c; begin hx = 0;
      for (i = 0; i < 4; i = i + 1) begin c = s[31 - 8 * i -: 8]; hx = (hx << 4) | ((c >= "A") ? (c - 55) : (c - 48)); end
    end
  endfunction
  localparam [15:0] I0 = hx(INIT0); localparam [15:0] I1 = hx(INIT1);
  wire inj = (INJECT == "YES"); wire z0, z30, z1, z31;
  lut4n u0(I0, A0, B0, C0, D0, z0, z30); lut4n u1(I1, A1, B1, C1, D1, z1, z31);
  wire cin = (CIN === 1'bz) ? 1'b0 : CIN;
  assign S0 = z0 ^ (cin & ~inj); wire c0 = z0 ? cin : (z30 & ~inj);
  assign S1 = z1 ^ (c0 & ~inj); assign COUT = z1 ? c0 : (z31 & ~inj);
endmodule
