// reference: two separate signed multipliers
module unpack_ref (input [7:0] a, input [7:0] b, input [11:0] c, output [18:0] pa, output [18:0] pb);
  wire [11:0] ce = (c == 12'h800) ? 12'h801 : c;
  wire signed [19:0] xa = $signed(a) * $signed(ce);
  wire signed [19:0] xb = $signed(b) * $signed(ce);
  assign pa = xa[18:0]; assign pb = xb[18:0];
endmodule
