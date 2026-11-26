function [S,ts]=mimoctf2dstep(num,den,tau,Ts,N)
%
% mimoctf2dstep     Computes the discrete-time step response
%                   of a continuous-time MIMO transfer function.
%
%                   Let a MIMO system (p outputs, m inputs) be specified by a continuous-time
%                   transfer function model
%
%                       Y(s) = G(s)*U(s)
%
%                       Y(s) = [y1(s) y2(s) ... yp(s)]'
%                       U(s) = [u1(s) u2(s) ... um(s)]'
%
%                              | g11(s) ... g1m(s) |
%                       G(s) = |   .         .     |
%                              | gp1(s) ... gpm(s) |
%
%                     The individual transfer functions are of the form
%
%                       yi(s) = gij(s)*uj(s) 
%                             = hij(s)*exp(-tau_ij*s)*uj(s)
%                             = hij(s)*uj(s-tau_ij)
%
%                       hij(s) = num_ij(s)/den_ij(s)
%
%                       den_ij(s) = a_ij(0)*s^na_ij + ... + a_ij(na_ij-1)*s + a_ij(na_ij)
%                       num_ij(s) = b_ij(0)*s^nb_ij + ... + b_ij(nb_ij-1)*s + b_ij(nb_ij)
%
%                     The time delay, tau_ij, must be non-negative and the transfer
%                     function can only be realized if the denominator degree is
%                     at least the same size as the numerator degree.
%
%                     The MIMO discrete-time step response coefficients are computed
%                     by computing the discrete-time step response
%                     coefficients for each input-output combination
%                     individually assuming a zero-order-hold of the
%                     inputs. This method is efficient as SISO systems with
%                     delays tend to have low dimension and are easily
%                     realized.
%   
% 
% Syntax: [S,ts]=mimoctf2dstep(num,den,tau,Ts,N)
%
%           num         :   p-times-m cell array with the numerator polynomials
%           den         :   p-times-m cell array with the denominator polynomials
%           tau         :   p-times-m matrix with the time delays
%           Ts          :   Sampling time period
%           N           :   Final sampling time (tfinal = N*Ts)
%
%           S           :   Step response matrices (p-times-m-times-(N+1))
%           ts          :   Corresponding time vector
%

[p,m]=size(tau);

S = zeros(p,m,N+1);
for i=1:m
    for j=1:p
        [s,ts]=sisoctf2dstep(cell2mat(num(j,i)),cell2mat(den(j,i)),tau(j,i),Ts,N);
        for k=1:N+1
            S(j,i,k)= s(k,1);
        end
    end
end
