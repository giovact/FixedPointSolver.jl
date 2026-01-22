dir=__EXAMPLES
npointsInt=1000
bound=6.0
tol=1e-6
tolsc=1e-4
damp=0.9
niter=5000

julia -t auto scripts/explore_multipleinits.jl --Model HopfieldExample \
                                                --param_loop %alpha% 0.002 0.15 step 0.005  linear \
                                                --param_loop %beta% 0.02 1.2 step 0.005 reciprocal \
                                                --OPnames M Q \
                                                    --OP0 1.0 1.0 \
                                                    --OP0 0.0 0.99 \
                                                --IntMethod LegendreQuadrature --npointsInt $npointsInt --bound $bound --tol $tol --maxiter $niter --damp $damp --dir $dir \
                                                --plot true

exit
julia -t auto scripts/critical_line.jl --Model HopfieldExample \
                                                --OP M 1.0 --OP Q 1.0 \
                                                --param_loop %beta% 0.01 0.01 1.0 reciprocal \
                                                --pscan %alpha% --pstart 1e-5 --delta_p 0.01 --p_extrema 0.0 0.16 --protocol follow --tol_scancrit $tolsc --increase true \
                                                --critical_line_name MRSG --function isRetrieval \
                                                --IntMethod LegendreQuadrature --npointsInt $npointsInt --bound $bound --tol $tol --maxiter $niter --damp $damp --dir $dir \
                                                --saveeachscan false 

julia -t auto scripts/critical_line.jl --Model HopfieldExample \
                                                --OP M 1.0 --OP Q 1.0 \
                                                --param_loop %beta% 0.01 0.01 1.0 reciprocal \
                                                --pscan %alpha% --pstart 1e-5 --delta_p 0.01 --p_extrema 0.0 0.1 --protocol follow --tol_scancrit $tolsc --increase true \
                                                --critical_line_name RMR --function isRetrievalStateDominant \
                                                --IntMethod LegendreQuadrature --npointsInt $npointsInt --bound $bound --tol $tol --maxiter $niter --damp $damp --dir $dir \
                                                --saveeachscan false
