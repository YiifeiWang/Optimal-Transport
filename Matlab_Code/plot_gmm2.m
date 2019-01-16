close all

opts=[];
opts.data=4;
gmmopt = [];
gmmopt.mu1    = [0.2;0.5]; 
gmmopt.mu2    = [0.6;0.7]; 
gmmopt.sigma1 = reshape([0.05, 0.03].^2, 1, 1, 2); 
gmmopt.sigma2 = reshape([0.03, 0.05].^2, 1, 1, 2); 
gmmopt.rho1   = [0.6,0.4];
gmmopt.rho2   = [0.3,0.7];
opts.gmmopt=gmmopt;
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

max_iters=[1e1,1e2,1e3,4e3];
epsilon=[1e-2,1e-3,1e-4,1e-5];

% call mosek
out1 = LP_mosek(model); 
pi1 = full(out1.pi);
for s=1:4
	% call sinkhorn
	algopt = [];
	algopt.add_coup = 1;
	algopt.epsilon = epsilon(s);
	algopt.maxiter = max_iters(s);
	out2 = LPER_shns(model,algopt);
	pi2 = out2.pi;

	subplot_tight(2,4,s); imshow(addframe(imfuse(-pi2, -pi1)));

	algopt = [];
	algopt.add_coup = 1;
	algopt.epsilon = epsilon(s);
	algopt.epsilon0 = epsilon(s)*10^5;
	algopt.maxiter = floor(max_iters(s)/5);
	out3 = LPER_shnsc(model,algopt);
	pi3 = out3.pi;

	subplot_tight(2,4,s+4); imshow(addframe(imfuse(-pi3, -pi1)));
	h=title(['$\epsilon$','=',mat2str(epsilon(s))]);
	set(h,'Interpreter','latex');
end
