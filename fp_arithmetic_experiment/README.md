# Experimental Pure C Arithmetic

This small example contains a testing and benchmarking function for an alternative implementation of Fp arithmetic for the prime p = `5*2**248 - 1`.
I do not use any assembly or specific instructions, other than those which the compiler finds.

The API is slightly different (I work with pointers more often than the original code) and the mathematics involved for the arithmetic follows some ideas I had with my collegue Thomas Pornin during Real World Crypto. It should be trivial to integrate this.

This code (aside from inversion and legendre) follows the Rust code. The timing seems a little better.

## Performance Gain

At the moment, it looks as though we get a roughly 3x speed up compared to the assembly when using my machine. I think there are some dumb choices I have made so over the next week I will try and improve this further.

## Building

To build the binary, I hope running `make` is enough. It works for my Mac, but maybe something different happens for ARM macs or linux machines... Forgive me, but cross platform building is not something I understand so this really is a "it works on my machine".

## Testing

I followed the SQIsign testing, so you run:

```
Jack: fp_arithmetic_experiment % ./fp_testing test 100000

--------------------------------------------------------------------------------------------------------

Testing field arithmetic over GF(p): 

  GF(p) equality tests ............................................ PASSED
  GF(p) addition tests ............................................ PASSED
  GF(p) subtraction tests ......................................... PASSED
  GF(p) montgomery conversion tests ............................... PASSED
  GF(p) multiplication tests ...................................... PASSED
  GF(p) squaring tests............................................. PASSED
  GF(p) inversion tests............................................ PASSED
  Square root, square tests........................................ PASSED
  ```

## Benchmarks

The following benchmarks are found using a 2.6 GHz 6-Core Intel Core i7 with turbo boost disabled

```
Jack: fp_arithmetic_experiment % ./fp_testing bench 10000

--------------------------------------------------------------------------------------------------------

Benchmarking GF(p) field arithmetic: 

  GF(p) addition runs in .......................................... 5 cycles, (5901499188998290977 ignore me)
  GF(p) subtraction runs in ....................................... 11 cycles, (12179369125404802345 ignore me)
  GF(p) multiplication runs in .................................... 24 cycles, (3586068188944577399 ignore me)
  GF(p) squaring runs in .......................................... 22 cycles, (3876598179343484886 ignore me)
  GF(p) inversion runs in ......................................... 5902 cycles, (17670238964497910354 ignore me)
  GF(p) sqrt runs in .............................................. 5624 cycles, (5452856236907160256 ignore me)
  Square checking runs in ......................................... 5890 cycles, (8869800425310996852 ignore me)
```

As a comparison, the following timings from the current repo

### Assembly

```
Benchmarking GF(p) field arithmetic: 

  GF(p) addition runs in .......................................... 13 cycles, (3024649066891134516 ignore me)
  GF(p) subtraction runs in ....................................... 11 cycles, (13394949140676670723 ignore me)
  GF(p) multiplication runs in .................................... 56 cycles, (8963423273615996848 ignore me)
  GF(p) inversion runs in ......................................... 14686 cycles, (5077583204304503167 ignore me)
  GF(p) sqrt runs in .............................................. 13920 cycles, (15383613233214440179 ignore me)
  Square checking runs in ......................................... 14660 cycles, (4391891120532673611 ignore me)
```

### Reference

```
Benchmarking GF(p) field arithmetic: 

  GF(p) addition runs in .......................................... 28 cycles, (4347325391616205458 ignore me)
  GF(p) subtraction runs in ....................................... 24 cycles, (4218286412216037907 ignore me)
  GF(p) multiplication runs in .................................... 107 cycles, (4231889596301059091 ignore me)
  GF(p) inversion runs in ......................................... 25655 cycles, (5007120656238929533 ignore me)
  GF(p) sqrt runs in .............................................. 24139 cycles, (14155235593077854046 ignore me)
  Square checking runs in ......................................... 25521 cycles, (8046875758828087190 ignore me)
```