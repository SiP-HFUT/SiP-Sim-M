**Note**: for *Loop_with_grating_v1.m* code, when simulating common Bragg gratings, using `x2 = [a2;r*a2+t*d2;t*a2+r*d2;d2];` together with `r = 1j * r` works fine.

But, when simulating a phase shifted Bragg grating, `x2 = [a2;r*a2+t*d2;t*a2-r*d2;d2];` and `r = -1j * r` works fine. This combination, however, doesn't work for common Bragg gratings.

The reason is unknown.