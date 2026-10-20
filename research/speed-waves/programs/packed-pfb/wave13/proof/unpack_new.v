// ONE DSP48E2, two 8-bit samples packed on A, one 12-bit coefficient on B; combinational (all registers off).
// The coefficient ROM never holds -2048 (symmetric range): c = -2048 is mapped to -2047 here and in the reference.
module unpack (input [7:0] a, input [7:0] b, input [11:0] c, output [18:0] pa, output [18:0] pb);
  wire [11:0] ce = (c == 12'h800) ? 12'h801 : c;
  wire [12:0] nc = -{ce[11], ce};
  wire [47:0] p;
  DSP48E2 #(.AREG(0), .BREG(0), .CREG(0), .MREG(0), .PREG(0), .ACASCREG(0), .BCASCREG(0), .OPMODEREG(0), .ALUMODEREG(0),
            .INMODEREG(0), .CARRYINREG(0), .CARRYINSELREG(0), .DREG(0), .ADREG(0)) u (
    .CLK(1'b0), .A({{3{a[7]}}, a, ~b[7], {10{b[7]}}, b}), .B({{6{ce[11]}}, ce}), .C({{17{nc[12]}}, nc, 18'd0}), .D(27'd0),
    .PCIN(48'd0), .ACIN(30'd0), .BCIN(18'd0), .OPMODE(9'b000110101), .ALUMODE(4'b0000), .INMODE(5'b00000),
    .CARRYINSEL(3'b000), .CARRYIN(1'b0), .CARRYCASCIN(1'b0), .MULTSIGNIN(1'b0),
    .CEA1(1'b1), .CEA2(1'b1), .CEAD(1'b1), .CEALUMODE(1'b1), .CEB1(1'b1), .CEB2(1'b1), .CEC(1'b1), .CECARRYIN(1'b1),
    .CECTRL(1'b1), .CED(1'b1), .CEINMODE(1'b1), .CEM(1'b1), .CEP(1'b1),
    .RSTA(1'b0), .RSTALLCARRYIN(1'b0), .RSTALUMODE(1'b0), .RSTB(1'b0), .RSTC(1'b0), .RSTCTRL(1'b0), .RSTD(1'b0),
    .RSTINMODE(1'b0), .RSTM(1'b0), .RSTP(1'b0), .P(p), .PCOUT());
  assign pb = p[18:0];
  assign pa = p[37:19] + p[18];
endmodule
