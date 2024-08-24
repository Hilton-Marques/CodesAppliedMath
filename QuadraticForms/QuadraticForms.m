classdef QuadraticForms < Scene
    properties
        m_A
    end
    methods
        function this = QuadraticForms(A)
            this = this@Scene(filename="q_forms", type="R2");
            this.m_A = A;
        end
        %Here we calculate C, such that A2 = C * A * C'
        function C = CalculateCongruentTransformation(this, A2)
            C = chol(A2)'/(chol(this.m_A)');
        end
        function Transform(this, A2)
            %define circle
            n = 100;
            theta = linspace(0,2*pi,n);
            x = [cos(theta); sin(theta)];

            C = this.CalculateCongruentTransformation(A2);
            C_1 = chol(A2)';
            C_2 = inv((chol(this.m_A)'));
            n = 50;
            for i = 1:n
                h = i/n;
                C_i_1 = this.Interp3DMatrix(C_1,h);
                C_i_2 = this.Interp3DMatrix(C_2,h);
                C_i = C_i_1 * C_i_2;
                A_i = C_i * this.m_A * C_i'; 
                T_i = inv(chol(A_i)');
                color_i = (this.m_blue + h * (this.m_red - this.m_blue));
                y = T_i * x;
                this.FillPolygon(y, color_i);
                this.get()
            end
            T_i = inv(chol(A2)');
            y = T_i * x;
            this.FillPolygon(y, color_i);
            this.save();
            this.exportFrame();
				end

				function A_i = Translate(this, A2,trans)
            %define circle
            n = 100;
            theta = linspace(0,2*pi,n);
            x = [cos(theta); sin(theta)];

            C = this.CalculateCongruentTransformation(A2);
            C_1 = chol(A2)';
            C_2 = inv((chol(this.m_A)'));
            n = 20;
						fac = 20;
						%trans = [1;0];
						h = 1/(n);
            for i = 1:n+1
                h_i = (i-1)*h;
                C_i_1 = this.LinearInterp3DMatrix(C_1,h_i);
                C_i_2 = this.LinearInterp3DMatrix(C_2,h_i);
                C_i = C_i_1 * C_i_2;
                A_i = C_i * this.m_A * C_i';

                T_i = inv(chol(A_i)');
                color_i = (this.m_blue + h_i * (this.m_red - this.m_blue));
                y = T_i * x;
								y = y + [fac*h_i*trans(1);trans(2)];
                this.FillPolygon(y, color_i);
            end
            T_i = inv(chol(A2)');
% 						y = T_i * x;
% 						y = y + [fac*trans(1);trans(2)];
% 						this.FillPolygon(y, color_i);
%             this.exportFrame();
				end

				function B_i = TranslateRiemmannian(this, B,trans)
					%define circle
					n = 100;
					theta = linspace(0,2*pi,n);
					x = [cos(theta); sin(theta)];

					%define n steps
					n = 20;
					t = linspace(0,1, n+1);
					h = t(2);
					sqrA = sqrtm(this.m_A);
					ABA = sqrA\B/sqrA; %A^-1/2*B*A^-1/2
					ABAt = ABA^h;
					temp = eye(2);
					fac = 20;
					%trans = [1;0];
					h = 1/(n);
					for i = 1:n+1
						h_i = (i-1)*h;
						B_i = sqrA * temp * sqrA;
						temp = temp * ABAt;
						%plot
						T_i = inv(chol(B_i)');
						color_i = (this.m_blue + h_i * (this.m_red - this.m_blue));
						y = T_i * x;
						y = y + [fac*(h_i)*trans(1);trans(2)];
						this.FillPolygon(y, color_i);
					end
					T_i = inv(chol(B)');
% 					y = T_i * x;
% 					y = y + [fac*trans(1);trans(2)];
% 					this.FillPolygon(y, color_i);
				end

				function B_i = TranslateEuclidean(this, B,trans)
					%define circle
					n = 100;
					theta = linspace(0,2*pi,n);
					x = [cos(theta); sin(theta)];

					%define n steps
					n = 20;
					t = linspace(0,1, n+1);
					h = t(2);
					ABAt = h*(B-this.m_A);
					temp = this.m_A;
					fac = 20;
					%trans = [1;0];
					h = 1/(n);
					for i = 1:n+1
						h_i = (i-1)*h;
						B_i = temp;
						temp = temp + ABAt;
						%plot
						T_i = inv(chol(B_i)');
						color_i = (this.m_blue + h_i * (this.m_red - this.m_blue));
						y = T_i * x;
						y = y + [fac*(h_i)*trans(1);trans(2)];
						this.FillPolygon(y, color_i);
					end
					T_i = inv(chol(B)');
% 					y = T_i * x;
% 					y = y + [fac*trans(1);trans(2)];
% 					this.FillPolygon(y, color_i);
				end

        function Show(this)
            L = chol(this.m_A)';
            T = inv(L);
            n = 100;
            theta = linspace(0,2*pi,n);
            x = [cos(theta); sin(theta)];
            y = T * x;
            this.FillPolygon(y, this.m_red)
        end
    end

end
