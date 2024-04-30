; -------------------------------------------------------------------------------------
; Test Amoeba a Simplex locally-convergent optimizer

pro test_amoeba

ftol=1e-4
initial=10.*randomu(seed,2)-5.0
print,'initial values =',initial
print
print,'The standard Nelder-Meade Simplex method (called here amoeba) can fail.'
print,'It suceeds for the Rosenbrock (without local minima) function.'
print
result=amoeba(ftol,function_name='himmelblau',p0=initial,scale=[5.,5.],ncalls=ncalls)
print,'Amoeba Himmelblau function'
print,'Global minimum at (3.,2.)=0.'
print,'parameters =',result
print,'function =',call_function('himmelblau',result)
print,'number of function calls=',ncalls
print

result=amoeba(ftol,function_name='rosenbrock',p0=initial,scale=[5.,8.],ncalls=ncalls)
print,'Amoeba Rosenbrock function'
print,'Global minimum at (1.,1.)=0.'
print,'parameters =',result
print,'function =',call_function('rosenbrock',result)
print,'number of function calls=',ncalls
print

ftol=1e-6
result=amoeba(ftol,function_name='rastrigin',p0=initial,scale=[5.,5.],ncalls=ncalls)
print,'Amoeba Rastrigin function'
print,'Global Optimum: -2, (0,0)'
print,'parameters =',result
print,'function =',call_function('rastrigin',result)
print,'number of function calls=',ncalls

end

; ----------------------------------------------------------------------------------
pro test_particle_swarm

population = 50
ftol=1e-2
print,'Rosenbrock function.  Global Optimum: 0,at (1,1)'  
particle_swarm,ftol,'rosenbrock',[-5.,-2.],[5.,8.],success,xmin,fxmin,fevals,PopSize=population,/print_results

ftol=1e-4
print,'Himmelblau function. Global Himmelblau Optimum: 0 (3,2)'
particle_swarm,ftol,'himmelblau',[-5.,-5.],[5.,5.],success,xmin,fxmin,fevals,PopSize=population,deflation=[[-3.78,-3.28],[-2.81,3.13],[3.58,-1.85]],/print_results

; True Optima: -1.031628453489877, (-0.08983,0.7126),
; (0.08983,-0.7126)
ftol = 1e-4
print,'Camel back'
particle_swarm,ftol,'camelback',[-5.,-5.],[5.,5.],success,xmin,fxmin,fevals,PopSize=population,/print_results

ftol=1e-4
print,'Goldstei-Price function. Global minimum (0.,1.)=3.0d0'
particle_swarm,ftol,'goldstein_price',[-2.,-2.],[2.,2.],success,xmin,fxmin,fevals,PopSize=population,MaxIt=100,/print_results

print,'Rastrigin function. Global Optimum: -2, (0,0)'
particle_swarm,ftol,'rastrigin',[-1.,-1.],[1.,1.],success,xmin,fxmin,fevals,PopSize=population,/print_results

print,'Griewank 2d. minimum at origin =0'
ftol = 1e-3
particle_swarm,ftol,'griewank_2d',[-600.,-600.],[600.,600.],success,xmin,fxmin,fevals,PopSize=population,/print_results

end

; ----------------------------------------------------------------------------
pro multi_minima_ps

wait=1
nminima = 4
ftol = 1e-4
; ----- 1st minimum (optimum)
particle_swarm,ftol,'himmelblau',[-5.,-5.],[5.,5.],success,xmin,fxmin,fevals,/parameter_plot
deflation = [xmin]

; -----
print,'Results'
print,'Parameters:',xmin
print,'Evaluation function value =',fxmin
print,'Number of evaluations =',fevals
; -----
kk=read_key(wait)  
iter = 0

if success eq 1 then begin
    while (success eq 1 and iter lt 10) do begin
        particle_swarm,ftol,'himmelblau',[-5.,-5.],[5.,5.],success $
                      ,xmin,fxmin,fevals,deflation=deflation,/parameter_plot
        if success eq 1 then begin
            print,'Minimum found'
        endif else begin
            print,'No minimum found'
        endelse
        print,'Results'
        print,'Parameters:',xmin
        print,'Evaluation function value =',fxmin
        print,'Number of evaluations =',fevals
        ; -----
        deflation = [deflation,xmin]
        kk=read_key(wait)
        iter = iter + 1
    endwhile    
endif

end

; ----------------------------------------------------------------------------
; Test functions

function himmelblau,x

; This is the Himmelblau Function
; Bound: X1=[-5,5], X2=[-5,5]
; Global Optimun: 0 (3,2)

f1 = (x(0)*x(0)+x(1)-11.0d0)
f2 = (x(0)+x(1)*x(1)-7.0d0)

f  = f1 + f2

return,f1*f1+f2*f2

end

function rosenbrock,x

;  The Rosenbrock Function is unimodal
;  Bound: X1=[-5,5], X2=[-2,8] 
;  Global Optimum: 0,at (1,1) is not easy to find because it is situated in a
;  valley with a flat bottom.

x1=x(0)
x2=x(1)
a = 100.0d0
f=a*(x2-x1*x1)^2+(1.0-x1)^2

return,f
 
end

function camelback,x

;  This is the Six-hump Camelback Function.
;  Bound: X1=[-5,5], X2=[-5,5]
;  True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
;
x1 = x(0)
x2 = x(1)
f = (4.0d0-2.1d0*x1^2+x1^(4.0d0/3.0d0))*x1^2+x1*x2+(-4.0d0+4.0d0*x2^2)*x2^2
    
return,f

end

function goldstein_price,x

; This is the Goldstein-Price Function
; Bound X1=[-2,2], X2=[-2,2]
; Global Optimum: 3.0,(0.0,-1.0)
;      
x1 = x(0)
x2 = x(1)
u1 = (x1 + x2 + 1.0d0)^2
u2 = 19.0d0 - 14.0d0*x1 + 3.0d0*x1^2 - 14.0d0*x2 + 6.0d0*x1*x2 +3.0d0*x2^2
u3 = (2.0d0*x1 - 3.0d0*x2)^2
u4 = 18.0d0 - 32.0d0*x1 + 12.0d0*x1^2 + 48.0d0*x2 -36.0d0*x1*x2 + 27.0d0*x2^2
u5 = u1 * u2
u6 = u3 * u4
f = (1.0d0 + u5) * (30.0d0 + u6)

return, f

end

function rastrigin,x

;
;  The Rastrigin Function has many hills and valleys
;  Bound: X1=[-1,1], X2=[-1,1]
;  Global Optimum: -2, (0,0)
 
x1 = x(0)
x2 = x(1)
f = x1^2 + x2^2 - cos(18.0d0*x1) - cos(18.0d0*x2)

return,f

end

function griewank_2d,x

;  This is the Griewank Function (2-D or 10-D)
;  Bound: X(i)=[-600,600], for i=1,2,...,10
;  Global Optimum: 0, at origin

d = 200.0

u1 = 0.0
u2 = 1.0
for j =0,1 do begin
    u1 = u1 + x(j)*x(j)/d
    u2 = u2 * cos(x(j)/sqrt(float(j+1)))
endfor
f = u1 - u2 + 1.

return,f

end

function griewank_10d,x

;  This is the Griewank Function (2-D or 10-D)
;  Bound: X(i)=[-600,600], for i=1,2,...,10
;  Global Optimum: 0, at origin

d = 4000.0

u1 = 0.0
u2 = 1.0
for j =0,9 do begin
    u1 = u1 + x(j)*x(j)/d
    u2 = u2 * cos(x(j)/sqrt(float(j+1)))
endfor
f = u1 - u2 + 1.

return,f

end

; ---------------------------------------------------------------------------
pro particle_swarm,ErrGoal,func,lwbnd,upbnd,success,xmin,fxmin,fevals,$
                   PopSize=PopSize,MaxIt=MaxIt,deflation=deflation, $
                   parameter_plot=parameter_plot,verbose=verbose,no_exceed=no_exceed, $
                   print_results=print_results

; Calling sequence
;    particle_swarm,1e-3,'func',[-2.,-2.],[2.,2.],xmin,fxmin,fevals $
;    ,PopSize = size of the swarm (default 20) $
;    ,MaxIt   = maximum number of iterations (default 1000) $
;    ,deflation = array of the position of parameters to be avoided $
;    ,parameter_plot = 1 if you want to plot the parameter plane (2
;                      parameters) $
;    ,verbose = 1 to print messages during the optimization
;    .print_results = 1 to print the final results  
; 
; Input
;    ErrGoal = desired accuracy (float or double)
;    func    = name of the function (string)
;    upbnd   = upper bound (float or double), array
;    lwbnd   = lower bound (float or double), array
;  
; Output
;    success = 0 (unsuccessful), 1 (successful)
;    xmin    = best parameters (array)
;    fxmin   = evaluation (merit, fitness) at the best position
;    fevals  = number of function evaluations
;
; Algorithm
;    Particle Swarm Optimization, dynamic inertia reduction variant
;
;    In the event that a particle falls outside the allowed domain, it is
;    given a new position randomly chosen inside the parameter domain
;
; Reference
;
;  Kennedy, J.; Eberhart, R. (1995). "Particle Swarm Optimization".
;  Proceedings of IEEE International Conference on Neural Networks.
;  IV. pp. 1942–1948.
;
;  Clerc, M. (2012). "Standard Particle Swarm Optimisation" (PDF).
;  HAL open access archive.  
;  https://hal.archives-ouvertes.fr/file/index/docid/764996
;  /filename/SPSO_descriptions.pdf
;
; Author
;    Wing-Fai Thi
;
; Version 1 22/3/2018 based on an old version from 2007
;
; ------------------------------
; --- Initializing variables

success  = 0    ; success flag
if not keyword_set(PopSize) then PopSize  = 20   ; size of the swarm
if not keyword_set(MaxIt) then MaxIt    = 1000  ; maximum number of iterations
iter     = 0    ; iterations' counter
fevals   = 0    ; function evaluations' counter
;maxw     = 1.2  ; maximum inertia weight's value
maxw     = 1.0  ; maximum inertia weight's value
minw     = 0.1  ; minimum inertia weight's value
weveryit = floor(0.75*MaxIt) ; inertia deer, step
c1       = 1.0  ; PSO parameter c1 Kennedy & Eberhart c1=c2=2.0
c2       = 1.0  ; PSO parameter c2
minc1    = 0.2
minc2    = 0.2 
inertdec = (maxw-minw)/float(weveryit) ; inertia weight's decrement
w        = maxw ; initial inertia weight
nparam      = n_elements(upbnd) ; number of parameters of the problem
gnrng    = dblarr(nparam) ; geometric mean of the parameters
epsilon  = 1d-30

if keyword_set(verbose) then begin
    print,'Initializing the variables'
    print,'Population Size =',Popsize
    print,'Maximum number of iterations =',MaxIt
endif

; ---

h     = 5
alpha = 0.9
wait=1
nexceed = 0
; --- Initializing swarm and velocity
popul = fltarr(nparam,Popsize)   ; nparam = column Popsize raw
bound = upbnd-lwbnd
for i=0,nparam-1 do popul(i,*) = randomu(seed,Popsize)*bound(i)+lwbnd(i)
vel   = randomu(seed,nparam,Popsize)

; --- Evaluate initial population
fpopul = fltarr(Popsize)
if keyword_set(verbose) then begin
    print,'Evaluate initial population'
endif
for i=0,Popsize-1 do begin
    fpopul(i)=call_function(func,popul(*,i))
    fevals = fevals + 1
endfor

; --- Initializing Best positions' matrix and the 
;     corresponding function values

bestpos  = popul  ; personal best
fbestpos = fpopul

; --- Finding the best particle in initial population

fbestpart = min (fpopul,g)
lastbpf   = fbestpart

; ---
window,0,retain=2
fbest = [1e10]
mean_fitness_past = 1e10

; --- Swarm evolution loop (start)
while (success eq 0) and (iter lt MaxIt) do begin

iter = iter + 1

; --- Update the value of the inertial weight w

;if iter le weveryit then w = maxw - (iter-1)*inertdec

; --- 
if (iter gt h) then begin
 if fbestpart ge fbest(iter-h) then begin
     w = alpha  * w
     c1 = alpha * c1
     c2 = alpha * c2
 endif    
endif

if (float(nexceed)/PopSize gt 0.2) then begin
    w = alpha * w
    c1 = alpha * c1
    c2 = alpha * c2
endif

w  = max(w,minw)
c1 = max(c1,minc1)
c2 = max(c2,minc2)


; --- Velocity update

if keyword_set(no_exceed) then begin
populvel = fltarr(nparam)
nexceed = 0
for i=0,Popsize-1 do begin
      ; repeat the velocity update until the new parameters stay into boundaries
      nrep = 0
      repeat begin
         nrep = nrep + 1 
         r1  = randomu(seed,nparam)
         r2  = randomu(seed,nparam)
         new_vel =  (w*vel(*,i) + c1*r1*(bestpos(*,i)-popul(*,i)) $
                               +  c2*r2*(bestpos(*,g)-popul(*,i)))/float(nrep) 
         populvel = popul(*,i)+new_vel
      endrep until (min(populvel-lwbnd) gt 0.0 and min(upbnd-populvel) gt 0.0) or nrep ge 5
    popul(*,i) = populvel
endfor
endif else begin

; --- Velocity update
a = fltarr(nparam,Popsize)
for i=0,PopSize-1 do a(*,i) = bestpos(*,g)
r1  = randomu(seed,nparam,PopSize)
r2  = randomu(seed,nparam,PopSize)
vel = w*vel + c1*R1*(bestpos-popul) + c2*R2*(a-popul)

; --- Swarm update
  
popul = popul + vel

endelse

; --- Assign a new random position and velocity if the particle lies
;     outside the boundaries
nexceed = 0
for i=0,PopSize-1 do begin
    if min(popul(*,i)-lwbnd) lt 0.0 or min(upbnd-popul(*,i)) lt 0.0 then begin
        popul(*,i) = randomu(seed,nparam)*bound+lwbnd
        vel(*,i)   = randomu(seed,nparam)
        nexceed = nexceed + 1
    endif
endfor   

; --- Evaluate the new swarm

for i=0,Popsize-1 do begin
    fpopul(i)=call_function(func,popul(*,i))
    if keyword_set(deflation) then begin
        dim_deflation = size(deflation)
        for j=0,dim_deflation(0)-1 do begin
            fpopul(i)=(fpopul(i)+epsilon)/(total((popul(*,i)-deflation(*,j))^2))
        endfor    
    endif    
    fevals = fevals + 1
endfor  
 
; --- Compute mean fitness

mean_fitness = median(fpopul)

; --- Updating the best position for each particle

WhereChanges = where(fpopul lt fbestpos,nchanges)
if nchanges gt 0 then begin
    fbestpos(WhereChanges)  = fpopul(WhereChanges)
    bestpos(*,WhereChanges) = popul(*,WhereChanges)
endif

; --- Updating index g (absolute best)

fbestpart = min(fbestpos,g)
fbest     = [fbest,fbestpart]

; --- show plot (in 2 dimension case)

if keyword_set(parameter_plot) then begin
   plot,popul(0,*),popul(1,*),psym=1,xr=[lwbnd(0),upbnd(0)],$
        yr=[lwbnd(1),upbnd(1)],/xs,/ys
   oplot,[bestpos(0,g)],[bestpos(1,g)],psym=4
;   kk=read_key(wait)
endif


if keyword_set(verbose) then begin
 print,'iter =',iter,' mean/best fitness =',mean_fitness,fbestpart,' w=',w,' exceed=',nexceed
 print,popul(*,g)
endif

; ----- Computes the normalized geometric range of the parameters

for i=0,nparam-1 do gnrng(i)=exp(mean(alog((max(popul(i,*))-min(popul(i,*)))/bound(i))))

; --- Checking stopping criterion

criterion = abs((mean_fitness-mean_fitness_past)/mean_fitness_past)

if criterion le ErrGoal then begin
   success = 1
endif else begin
   lastbpf = fbestpart
   mean_fitness_past = mean_fitness
endelse

endwhile


; --- Swarm evolution loop (end)

; --- Output arguments
xmin  = popul(*,g)
fxmin = fbestpos(g)

; -----
if keyword_set(print_results) then begin
   print,'Results'
   print,'Parameters:',xmin
   print,'Evaluation function value =',fxmin
   print,'Number of evaluations =',fevals
   print
endif

end
