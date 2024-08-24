classdef RobustEkf < Scene
	properties
		m_H
		m_b
		m_psi
		m_iter
		m_rho
		m_prior
		m_post
		m_h
		m_stdy
		m_std_x
		m_xc
	end

	methods
		function this = RobustEkf(options)
			arguments
				options.H = 1
				options.b = 1
				options.method = "l2"
				options.iters (1,1) double = 10
				options.h = @(x) 1
				options.stdx = 1
				options.x_c = 1
			end
			this = this@Scene();
			H = options.H;
			b = options.b;
			method = options.method;
			iters = options.iters;
			h = options.h;
			stdx = options.stdx;
			x_c = options.x_c;


			this.m_prior = this.CreateGaussian(x_c, stdx);
			this.m_xc = x_c;
			this.m_std_x = stdx;
			this.m_H = H;
			this.m_b = b;
			this.m_iter = iters;
			this.m_h = h;
			switch method
				case 'l2'
					this.m_psi = @(x) ones(size(x));
					this.m_rho = @(x) x.^2;
				case 'huber'
					this.m_rho = @(x) this.huberRho(x, 0.1);
					this.m_psi = @(x) this.huberPsi(x, 1);
				case 'tukey'
					this.m_psi = @(x) tukey(x, 1);
				case 'bisquare'
					this.m_psi = @(x) bisquare(x, 1);
				otherwise
					error('Unknown method');
			end
		end

		function [new_x, C] = CalculatePosterior(this, z, std_y)
			z_bar = this.m_h(this.m_xc);
			dz = z - z_bar;
			C = [[this.m_std_x*this.m_std_x, 0]; [0, std_y * std_y]];
			Cinv = inv(C);
			H = [1;this.m_H];
			b = [0;dz];
			P_i = H'*Cinv*H;
			new_x = P_i^-1*H'*Cinv*b;
			new_xc = new_x + this.m_xc;
		end

		function [H, b] = GetCanonicalSpace(this, z, std_y)
			z_bar = this.m_h(this.m_xc);
			C = [[this.m_std_x*this.m_std_x, 0]; [0, std_y * std_y]];
			L  = chol(C);
			dz = z - z_bar;			
			H = [1; this.m_H];
			b = [0;dz];
			H = inv(L')*H;
			Cinv = eye(2);
			b = inv(L')*b;
		end

		function [x, z] = CalculateRobustPosterior(this, z, std_y)
			[H, b] = this.GetCanonicalSpace(z,std_y);
			x = pinv(H)*b;
			obj = IRWLS(H, b, "tukey",2.2);
			x = obj.Solve();
			z = obj.Show(x);
		end

		function ShowDistributionsPosterior(this, z, std_y, show_transformed,show_robust)
			arguments
				this
				z
				std_y
				show_transformed = false
				show_robust = false
			end
			axis normal
			[new_x, C] = this.CalculatePosterior(z, std_y);
			[~, contour] = this.CalculateRobustPosterior(z, std_y);
			%writematrix(contour)
			z_bar = this.m_h(this.m_xc);
			dz = z - z_bar;
			Cinv = inv(C);
			H = [1; this.m_H];
			b = [0;dz];
			if (show_transformed == true)
				L  = chol(C);
				H = inv(L')*H;
				Cinv = eye(2);
				b = inv(L')*b;
				new_x = pinv(H)*b;
				Y = H*new_x;
				v = L' * Y;
				c = [this.m_xc; z_bar];
				this.arrow([c;0],[c + Y;0],'color', this.m_green,'tipWidth', 0.060,'stemWidth', 0.035);
				this.arrow([c;0],[c + v;0],'color', this.m_blue,'tipWidth', 0.060,'stemWidth', 0.035);
				axis equal
			end
			axis equal
			P_i = H'*Cinv*H;
			new_xc = new_x + this.m_xc;
			
			%% Show priors
			c = [this.m_xc; z];
			margin_left = 8;
			margin_right = 4;
			x = linspace(new_xc - margin_left, new_xc + margin_right, 300);
			this.ShowGaussian(this.m_xc, this.m_std_x, x=x, trans=-z,color="#699C52");
			this.ShowGaussian(new_xc, 1/sqrt(P_i), trans=z, x=x,color=this.m_blue);
			if show_robust 
				fill3(contour(1,:)+c(1), contour(2,:)+c(2),0*contour(1,:),'red','FaceAlpha',0.4);
			else
				this.ShowEllipse(Cinv,'red', c );
			end
			
			%% Show measurement function
			y = (H(2)/H(1))*(x - this.m_xc) + z_bar;
			plot(x,y,"color",'#1C758A','linewidth',1.5);
			%% Show lines
			yline(z, '--','LineWidth',1)
			xline(new_x + this.m_xc,'--','LineWidth',1)
			yline(z_bar,'--','LineWidth',1)
			xline(this.m_xc,'--','LineWidth',1)
			%% Important points
			c = [this.m_xc; z_bar];
			new_Y = H*P_i^-1*H'*Cinv*b;
			plot(new_xc, z , 'o','markersize',6,'MarkerFaceColor','black');
			plot(this.m_xc, z, 'o','markersize',6,'MarkerFaceColor','black');
			plot(this.m_xc, -z, 'o','markersize',6,'MarkerFaceColor','black');
			plot(new_xc, new_Y(2) + z_bar, 'o','markersize',6,'MarkerFaceColor','black');
			this.arrow([c; 0.4],[c + b; 0.4],'color', '#699C52','tipWidth', 0.060,'stemWidth', 0.035);
			this.arrow([c+b;0],[c + (new_Y);0],'color', '#644172','tipWidth', 0.060,'stemWidth', 0.035);
			this.arrow([c; 0.4],[c + (new_Y);0],'color', this.m_blue,'tipWidth', 0.060,'stemWidth', 0.035);
			xlabel('$x$','interpreter','latex','FontSize',22)
			ylabel('$y$','interpreter','latex','FontSize',22)
			ylim([-z,2*z]);
			camlight
			
			this.exportFrame(filename="ekf_normalized_robust");
		end

		function h = ShowGaussian(this, x_c, std_x, options)
			arguments
				this
				x_c
				std_x
				options.fac = 10;
				options.trans = 0;
				options.x = linspace(-1,1,30);
				options.color = 'black';
			end
			trans = options.trans;
			fac = options.fac;
			x = options.x;
			color = options.color;

			gaussian = @(x) (1/(std_x * sqrt(2*pi))) * ...
											exp(((x - x_c).^2 ./ (std_x^2)) * -0.5);			
			h = plot(x, fac*gaussian(x) + trans, ...
				'color', color,'LineWidth',1.5); %prior
		end

		function x = Solve(this)
			% Solve the problem
			S = pinv(this.m_H);
			x = S*this.m_b;
			% Compute the residuals
			r = this.m_b - this.m_H*x;
			for i = 1:this.m_iter
				% Compute the weights
				w = this.m_psi(r);
				% Compute the weighted least squares solution
				S = diag(w)*this.m_H;
				x = pinv(S)*w*this.m_b;
				% update the residuals
				r = this.m_b - this.m_H*x;
			end
		end

		function Show(this, x)
			% Create original function
			r = this.m_H*x - this.m_b;
			z_value = this.J(r);

			% Create Contour plot
			n = 100;
			x_inv = linspace(-1,1,n);
			y_inv = linspace(-1,1,n);
			[x_grid, y_grid] = meshgrid(x_inv,y_inv);
			x = [x_grid(:), y_grid(:)];
			values = this.J(x);
			z = reshape(values, n, n);
			contourf(x_grid, y_grid, z,'FaceAlpha',0.5,'FaceColor','blue');
		end
		function res = J(this, r)
			% Compute the objective function
			if size(r,2) == 1
				r = r';
			end
			n = size(r,1);
			res = zeros(n,1);
			for i = 1:n
				res(i) = sum(this.m_rho(r(i,:)));
			end
		end
	end

	methods (Static)
		function x = huberRho(x, c)
			% Huber function
			x = abs(x);
			x(x <= c) = 0.5*x (x <= c).^2;
			x(x > c) = c*(x(x > c) - 0.5*c);
		end

		function x = huberPsi(x, c)
			% Huber function
			x = abs(x);
			x(x <= c) = x(x <= c) ^ 2;
			x(x > c) = c;
		end

		function x = tukey(x, c)
			% Tukey function
			x = abs(x);
			x(x <= c) = (c^2/6)*(1 - (1 - (x(x <= c)/c).^2).^3);
			x(x > c) = c^2/6;
		end

		function x = bisquare(x, c)
			% Bisquare function
			x = abs(x);
			x(x <= c) = (1 - (x(x <= c)/c).^2).^2;
			x(x > c) = 0;
		end

		function gaussian = CreateGaussian(x_c, std)
			gaussian = @(x) (1/(std_x * sqrt(2*pi))) * ...
											exp(((x - x_c).^2 ./ (std^2)) * -0.5);
		end
	end
end