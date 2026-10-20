// as-4: one DSP48E1, all registers off: A = a*2^17 + b + 2^16 (25 bits used), B = c (11-bit two's complement, sign-
// extended), C = -c*2^16; P = (a*c)*2^17 + b*c. Extraction: yb = P[16:0] (signed), ya = P[33:17] + P[CB].
// CB = 16 is correct; CB = 15 is the deliberately broken copy.
module unpack17 #(parameter CB = 16) (input [7:0] a, input [7:0] b, input [10:0] c, output [16:0] ya, output [16:0] yb);
  wire [29:0] A = {{5{a[7]}}, a, ~b[7], {8{b[7]}}, b};
  wire [17:0] B = {{7{c[10]}}, c};
  wire signed [47:0] cneg = -($signed({{37{c[10]}}, c}) <<< 16);
  wire [47:0] P;
  DSP48E1 #(.AREG(0), .BREG(0), .ACASCREG(0), .BCASCREG(0), .CREG(0), .MREG(0), .PREG(0), .OPMODEREG(0), .ALUMODEREG(0),
            .INMODEREG(0), .CARRYINREG(0), .CARRYINSELREG(0)) u (
    .CLK(1'b0), .A(A), .B(B), .C(cneg), .D(25'd0), .PCIN(48'd0), .ACIN(30'd0), .BCIN(18'd0),
    .OPMODE(7'b0110101), .ALUMODE(4'b0000), .INMODE(5'b00000), .CARRYINSEL(3'b000), .CARRYIN(1'b0),
    .CARRYCASCIN(1'b0), .MULTSIGNIN(1'b0),
    .CEA1(1'b1), .CEA2(1'b1), .CEAD(1'b1), .CEALUMODE(1'b1), .CEB1(1'b1), .CEB2(1'b1), .CEC(1'b1), .CECARRYIN(1'b1),
    .CECTRL(1'b1), .CED(1'b1), .CEINMODE(1'b1), .CEM(1'b1), .CEP(1'b1),
    .RSTA(1'b0), .RSTALLCARRYIN(1'b0), .RSTALUMODE(1'b0), .RSTB(1'b0), .RSTC(1'b0), .RSTCTRL(1'b0), .RSTD(1'b0),
    .RSTINMODE(1'b0), .RSTM(1'b0), .RSTP(1'b0), .P(P), .PCOUT());
  assign yb = P[16:0];
  assign ya = P[33:17] + P[CB];
endmodule
module top (input [7:0] a, input [7:0] b, input [10:0] c, output [16:0] ya, output [16:0] yb, output [16:0] ya_br, output [16:0] yb_br);
  unpack17 #(16) ok (.a(a), .b(b), .c(c), .ya(ya), .yb(yb));
  unpack17 #(15) br (.a(a), .b(b), .c(c), .ya(ya_br), .yb(yb_br));
endmodule
