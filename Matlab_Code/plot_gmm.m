close all

opts=[];
opts.data=4;
model=model_unified(opts);
m = model.m;
cst = model.cst;
mu = cst(1:m);
nu = cst(m+1:2*m);

index = [0:1/(m-1):1]';
figure('Renderer', 'painters', 'Position', [10 10 900 300]);
plot(index,mu); 
h = title('$\mu$');
set(h,'Interpreter','latex');
hold on
plot(index,nu);
h = title('$\nu$');
set(h,'Interpreter','latex');

figure('Renderer', 'painters', 'Position', [10 10 900 300]);

max_iters=[2e3,2e3,2e3,2e3];
epsilon=[1,1e-1,1e-2,1e-3];

% call mosek
out1 = LP_mosek(model); 
pi1 = full(out1.pi);
for s=1:4
	% call sinkhorn
	algopt = [];
	algopt.add_coup = 1;
	algopt.epsilon = epsilon(s);
	algopt.maxiter = max_iters(s);
	out2 = LPER_sinkhorn(model,algopt);
	pi2 = out2.pi;

	subplot_tight(1,4,s); imshow(addframe(imfuse(-pi2, -pi1)));
	h=title(['$\epsilon$','=',mat2str(epsilon(s))]);
	set(h,'Interpreter','latex');
end
