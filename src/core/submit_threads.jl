function submit_threads(job::Job{Body}, n_chunks::Int)
    point_decomp = defaultdist(job.spatial_setup.n_points, n_chunks)
    body_chunks = chop_body_threads(job.spatial_setup, job.time_solver, point_decomp)
    return nothing
end
