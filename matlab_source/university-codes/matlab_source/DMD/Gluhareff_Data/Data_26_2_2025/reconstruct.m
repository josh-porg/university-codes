function reconstruction = reconstruct(modes,timeDynamics)
%RECONSTRUCT Reconstructs based on modes and time dynamics
%   To Be used in conjuction with DMD or POD
reconstruction = modes * timeDynamics;
end

