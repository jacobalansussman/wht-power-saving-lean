// full-width helper: the same lines as unpack_new.v, the DSP replaced by its arithmetic at 64 bits
module unpack_helper (input [7:0] a, input [7:0] b, input [11:0] c, output [18:0] pa, output [18:0] pb);
  wire [11:0] ce = (c == 12'h800) ? 12'h801 : c;
  wire signed [63:0] A = $signed({a, ~b[7], {10{b[7]}}, b});
  wire signed [63:0] p = A * $signed(ce) - ($signed(ce) <<< 18);
  assign pb = p[18:0];
  assign pa = p[37:19] + p[18];
endmodule
