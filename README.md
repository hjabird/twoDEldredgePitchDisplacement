OpenFOAM® and OpenCFD® are registered trademarks of OpenCFD Limited, the producer OpenFOAM software. All registered trademarks are property of their respective owners.
This offering is not approved or endorsed by OpenCFD Limited, the producer of the OpenFOAM software and owner of the OPENFOAM® and OpenCFD® trade marks.
http://openfoam.org/

Based on a smooth pitching motion function from Eldridge 2009.

In Julia the function looks like (copied from K Ramesh's UNSFlow):
```
immutable EldUpDef <: MotionDef
    amp :: Float64  // Amplitude
    K :: Float64    // pitch rate non-dimensional parameter
    a :: Float64    // relation trasition speed
end

function (eld::EldUpDef)(t)
    sm = pi*pi*eld.K/(2*(eld.amp)*(1 - eld.a))
    t1 = c / U // Modified to consider time
    t2 = t1 + ((eld.amp)/(2*eld.K)) * c / U
    ((eld.K/sm)*log(cosh(sm*(t - t1))/cosh(sm*(t - t2))))+(eld.amp/2)
end
```
