!-------------------------------------------------------------
!
!       include file for nucleon fsi history ( nucleonfsihist.h )
!
!       common block to store nucleon fsi history for writing to output file later
!
!-------------------------------------------------------------
        integer maxnucleonvert           ! maximum number of vertices
        parameter(maxnucleonvert=200)

        integer maxnucleonstep           ! maximum number of steps
        parameter(maxnucleonstep=10000)
  
!
!       NFnvert        : number of vertices
!
!       NFiflag(i)     : 4-digit flag for interaction type at i-th vertex, in the form "BNTP":
!                         N: charge nucleon propagated through nucleus (0 = neutron, 1 = proton)
!                         T: charge "target" nucleon the interaction is taking place on
!                         P: scattering process:
!                            P=0: start tracking of nucleon (i.e. gets "created")
!                            P=1: elastic scattering
!                            P=2: single pion production
!                            P=3: double pion production
!                            P=4: stop tracking of nucleon (i.e. leaves nucleus)
!                         B: Pauli blocking flag (0 = not blocked, 1 = interaction was Pauli blocked
!                            and actually did not take place)
!                         Examples:
!                          - 103 means double pion production when a proton scattered on a neutron
!                          - 1011 means elastic scattering of a neutron on a proton did not take
!                            place due to Pauli blocking
!                         For P=0 and P=4, "T" is without meaning and always set to 0.
!       NFx(i)         : x-component of i-th vertex position inside nucleus
!       NFy(i)         : y-component of i-th vertex position inside nucleus
!       NFz(i)         : z-component of i-th vertex position inside nucleus
!       NFpx(i)        : x-component of momentum of nucleon leaving the i-th vertex
!       NFpy(i)        : y-component of momentum of nucleon leaving the i-th vertex
!       NFpz(i)        : z-component of momentum of nucleon leaving the i-th vertex
!       NFe(i)         : energy of nucleon leaving the i-th vertex
!       NFfirststep(i) : first step index of this track (to obtain the CMS energies for each step)
!
!       NFnstep        : number of steps
!
!       NFecms2(k)     : CMS energy squared of collision at k-th step (i.e. before interacting).
!                        The sign of this value indicates the charge of the target nucleon:
!                         NFecms2 > 0: proton,  NFecms2 < 0: neutron (same as "T" in NFiflag)
!       NFptot(k)      : total probability at k-th step (for testing only, will be removed)
!
!       Remarks:
!        - a "vertex" is actually better described as a start, end or scattering point of a track
!        - at each scattering point, the first nucleon will be followed in the same track, whilst the
!          second one will create a new track
!        - each track consists of a series of consecutive vertices. The first vertex has P=0, the
!          last P=4. In between may be any number (including 0) vertices where an actual scattering
!          took place (P=1,2,3).
!        - it is not possible (and not needed) to connect the second track of a scattering vertex
!          with the original one. Note that "first" and "second" is purely arbitrary. For nucleon
!          FSI uncertainties, only the probabilities of the scattering processes have to be
!          calculated, so it is not important to know which tracks belong to each other.
!
!

        integer*4   NFnvert, NFiflag,NFnstep,NFfirststep,NReweightFlag,NFnucresflg,NFiflagStep,NFStepNumber,RW_step_i
        real*4      NFx,NFy,NFz,NFpx,NFpy,NFpz,NFe,NFecms2,NFiPel,NFiPsp,NFiPdp,
     &   NFptot,NFxStep,NFyStep,NFzStep,NFpxStep,NFpyStep,NFpzStep,NFeStep,NFfirststepStep,NFiRhon
	common /nucleonfsihist/ NFnvert, NFiflag(maxnucleonvert),NReweightFlag,
     &   NFnucresflg,NFiflagStep(maxnucleonstep),NFStepNumber(maxnucleonstep),RW_step_i,
     &   NFx(maxnucleonvert),NFy(maxnucleonvert),NFz(maxnucleonvert),
     &   NFpx(maxnucleonvert),NFpy(maxnucleonvert),NFpz(maxnucleonvert),
     &   NFe(maxnucleonvert),NFfirststep(maxnucleonvert),
     &   NFnstep, NFecms2(maxnucleonstep),NFiPel(maxnucleonstep),NFiPsp(maxnucleonstep),
     &   NFiPdp(maxnucleonstep),NFptot(maxnucleonstep),
     &   NFxStep(maxnucleonstep),NFyStep(maxnucleonstep),NFzStep(maxnucleonstep),
     &   NFpxStep(maxnucleonstep),NFpyStep(maxnucleonstep),NFpzStep(maxnucleonstep),
     &   NFeStep(maxnucleonstep),NFfirststepStep(maxnucleonstep),NFiRhon(maxnucleonstep)

