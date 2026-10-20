// ccu2c_model.v (as-5): the ECP5 CCU2C carry cell with the semantics of Yosys share/yosys/lattice/ccu2c_sim.vh
// (S = LUT4 ^ cin; cout = LUT4 ? cin : LUT2(A,B) = INIT[3:0]). LUTs are written as ?: mux trees (as Yosys' own
// LUT4 model), so an input the function does not depend on does not turn the output into X in simulation.
module lut4m(input [15:0] I, input A, B, C, D, output Z);
  wire [7:0] s3 = D ? I[15:8] : I[7:0]; wire [3:0] s2 = C ? s3[7:4] : s3[3:0];
  wire [1:0] s1 = B ? s2[3:2] : s2[1:0]; assign Z = A ? s1[1] : s1[0];
endmodule
module CCU2C(input CIN, input A0, B0, C0, D0, A1, B1, C1, D1, output S0, S1, output COUT);
  parameter [15:0] INIT0 = 16'h0000; parameter [15:0] INIT1 = 16'h0000;
  parameter INJECT1_0 = "YES"; parameter INJECT1_1 = "YES";
  wire l4_0, l2_0, l4_1, l2_1;
  lut4m u0(INIT0, A0, B0, C0, D0, l4_0); lut4m v0({12'd0, INIT0[3:0]}, A0, B0, 1'b0, 1'b0, l2_0);
  lut4m u1(INIT1, A1, B1, C1, D1, l4_1); lut4m v1({12'd0, INIT1[3:0]}, A1, B1, 1'b0, 1'b0, l2_1);
  wire gc0 = (INJECT1_0 == "YES") ? 1'b0 : CIN; assign S0 = l4_0 ^ gc0;
  wire g0 = (INJECT1_0 == "YES") ? 1'b0 : l2_0; wire cout_0 = l4_0 ? CIN : g0;
  wire gc1 = (INJECT1_1 == "YES") ? 1'b0 : cout_0; assign S1 = l4_1 ^ gc1;
  wire g1 = (INJECT1_1 == "YES") ? 1'b0 : l2_1; assign COUT = l4_1 ? cout_0 : g1;
endmodule
