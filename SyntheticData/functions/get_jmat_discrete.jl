### function for discretising continuous time Jacobian
function get_jmat_discrete(jmat_cont)
    eig = eigen(jmat_cont)
    P = eig.vectors
    v_eig = eig.values
    jmat_discrete = P*diagm(exp.(v_eig))*inv(P)
    jmat_discrete = real.(jmat_discrete)
    return jmat_discrete
end