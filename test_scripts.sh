
dir="__TEST"
nthreads=8

rm -rf $dir

julia -t $nthreads scripts/explore_multipleinits.jl --Model HopfieldExample \
                                                --param_loop %alpha% 0.04 0.15 step 0.03  linear \
                                                --param_loop %beta% 0.05 1.0 length 15 reciprocal \
                                                --OPnames M Q \
                                                    --OP0 1.0 1.0 \
                                                    --OP0 0.0 0.99 \
                                                --IntMethod LegendreQuadrature --npointsInt 1000 --bound 7.0 --tol 1e-6 --maxiter 2000 --damp 0.9 --dir $dir \
                                                --plot false


julia -t $nthreads scripts/critical_line.jl --Model HopfieldExample  \
                                                --OP M 1.0 --OP Q 1.0 \
                                                --param_loop %beta% 0.05 0.6 step 0.05 reciprocal \
                                                --pscan %alpha% --pstart 0.01 --delta_p 0.01 --p_extrema 0.0001 0.16 --protocol follow --tol_scancrit 1e-4 --increase true \
                                                --critical_line_name Retrieval --function isRetrieval \
                                                --IntMethod LegendreQuadrature --npointsInt 1000 --bound 7.0 --tol 1e-6 --maxiter 2000 --damp 0.9 --dir $dir \
                                                --saveeachscan false 


julia -t $nthreads scripts/parallel_scans.jl --Model HopfieldExample \
                                                --param_scan %alpha% 0.04 0.15 step 0.03  linear \
                                                --param_loop %beta% 0.05 1.0 length 15 reciprocal \
                                                --OPnames M Q \
                                                --OP0 1.0 0.7 \
                                                --IntMethod LegendreQuadrature --npointsInt 1000 --bound 7.0 --tol 1e-6 --maxiter 2000 --damp 0.9 --dir $dir \








exit






