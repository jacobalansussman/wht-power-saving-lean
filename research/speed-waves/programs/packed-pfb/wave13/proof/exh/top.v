module top (input [7:0] a, input [7:0] b, input [11:0] c, output [18:0] pa_n, output [18:0] pb_n,
            output [18:0] pa_r, output [18:0] pb_r, output [18:0] pa_x, output [18:0] pb_x);
  unpack        n (.a(a), .b(b), .c(c), .pa(pa_n), .pb(pb_n));
  unpack_ref    r (.a(a), .b(b), .c(c), .pa(pa_r), .pb(pb_r));
  broken_unpack x (.a(a), .b(b), .c(c), .pa(pa_x), .pb(pb_x));
endmodule
