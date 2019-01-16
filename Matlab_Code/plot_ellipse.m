opts.data = 2;
opts.size = 50;
model = model_unified(opts);
out = LP_gurobi(model);
model.x = out.pi;
vector_field(model)

