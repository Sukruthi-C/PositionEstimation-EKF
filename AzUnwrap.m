% AzUnwrap
%
% Inputs: Az, clearvar
% Outputs: AzOut
% Persistent Variables: kwrap, prevAz
%
% Description: Function that can be used to unwrap the azimuth angle
% measurements obtained from the estimated position components. kwrap
% represents the amount of times the angle has been wrapped (can be
% negative) and prevAz saves the last unwrapped azimuth value (note that
% the initial value of prevAz is the initial azimuth measurement).
%
% Usage Note 1: 
% Let Az be the azimuth angle in radians that you obtain from using
% the measurement nonlinear model with the current position estimate. Then
% AzOut = AzUnwrap(Az, 0) yields the unwrapped azimuth angle in radians
% AzOut.
%
% Usage Note 2:
% Before starting any estimation procedure (going inside a loop) use
% AzUnwrap(0, 1) to reset the persistent variables. Otherwise, they may
% keep the values from previous calls.
%
function AzOut = AzUnwrap(Az, clearVar)

    persistent kwrap
    persistent prevAz
    
    if clearVar == 1
       kwrap = [];
       prevAz = [];
       return
    end
    
    if isempty(kwrap)
        kwrap = 0;
        prevAz = -1.5708;
    end
    
    if ((Az + 2*kwrap*pi) - prevAz <= -pi)
       kwrap = kwrap + 1;
    elseif ((Az + 2*kwrap*pi) - prevAz >= pi)
       kwrap = kwrap - 1;
    end
    AzOut = Az + 2*kwrap*pi;
    prevAz = AzOut;

end