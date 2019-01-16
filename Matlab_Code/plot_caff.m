opts.data = 3;
opts.size = 100;
model = model_unified(opts);
out = LP_gurobi(model);
model.x = out.pi;
vector_field(model)

