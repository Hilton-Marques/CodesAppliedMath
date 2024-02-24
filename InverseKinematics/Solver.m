classdef Solver < Scene
    properties
        m_d1;
        m_d2;
        m_d3;
        m_target;
        m_alpha = 0.01;
        m_R1;
        m_R2;
        m_R3;
        m_handles = [];
        m_origin = [0;0;0.25]
    end
    methods
        function this = Solver(d1, d2, d3, target)
            this = this@Scene('inverse_kinematics.gif',"IK");
            this.m_d1 = d1;
            this.m_d2 = d2;
            this.m_d3 = d3;
            this.m_target = target;
            this.m_R1 = SO3(eye(3));
            this.m_R2 = SO3(eye(3));
            this.m_R3 = SO3(eye(3));
            view(38,19);
            axis tight
            this.show();
            this.exportFrame('first');
            this.get();
            %this.check_loss_grad();
            this.Compute();
        end
        function show(this)
            delete(this.m_handles);
            link1 = this.m_origin;
            link2 = this.m_R1.act(this.m_d1);
            link3 = this.m_R2.act(this.m_R1.act(this.m_d2));
            link4 = this.m_R3.act(this.m_R2.act(this.m_R1.act(this.m_d3)));
            links = [link1,link2,link3,link4];
            h = this.drawPlan(link1,0,this.m_blue);
            h = [this.drawLandmarks(this.m_target), h];
            h(end+1) = this.plotText(this.m_target, '$\mathbf{y}^*$');
            positions = zeros(3,4);
            pi = [0;0;0];
            for i = 1:4
                pj = pi + links(:,i);
                h = [this.drawArm(pi, pj), h];
                pi = pj;
                positions(:,i) = pj;
                if (i < 4)
                h = [this.drawSphere(0.155,pj,[]), h];
                end
            end
            pGarra = positions(:,end);
            h = [this.drawPlan(links(:,end),0,[.8,.8,.8],0.5,[0;0;0],pGarra), h];
            this.m_handles = h;
        end
        function y = forward(this,R1, R2, R3)
            y = R1.act(this.m_d1) + R2.act(R1.act(this.m_d2)) ...
               + R3.act(R2.act(R1.act(this.m_d3))) + this.m_origin;
        end
        function f = loss(this, R1, R2, R3)
            diff = this.forward(R1, R2, R3) - this.m_target;
            f = dot(diff, diff);
        end
        function setTarget(this, target)
            this.m_target = target;
        end
        function Compute(this)
            iter = 0;
            e = this.loss(this.m_R1, this.m_R2, this.m_R3);
            iters = 200;
            while (iter < iters)
                grad = -this.m_alpha * this.loss_grad(this.m_R1, ...
                                                      this.m_R2, ...
                                                      this.m_R3);
                t1 = grad(1:3)';
                t2 = grad(4:6)';
                t3 = grad(6:9)';
                this.m_R1 = this.m_R1 * t1;
                this.m_R2 = this.m_R2 * t2;
                this.m_R3 = this.m_R3 * t3;
                iter = iter + 1;
                if (iter < 10 || mod(iter,10) == 0)
                    this.show();
                    pause(0.1);
                    this.get();
                end
                e(end+1) = this.loss(this.m_R1, this.m_R2, this.m_R3);
                if e(end) > e(end-1)
                    break;
                end
            end
            y = this.forward(this.m_R1,this.m_R2,this.m_R3);
            norm(y);
            this.plotErrors(1:size(e,2), e);
        end
        function res = check_loss_grad(this)
            tau1 = SO3.getRandomVector()';
            tau2 = SO3.getRandomVector()';
            tau3 = SO3.getRandomVector()';
            R1 = SO3(SO3.getRandomElement());
            R2 = SO3(SO3.getRandomElement());
            R3 = SO3(SO3.getRandomElement());
            tau = [tau1;tau2;tau3];
            for i = 1:9
                y = this.loss(R1,R2,R3);
                diff = this.loss(R1 * tau1, R2 * tau2, R3 * tau3) - y;
                d = this.loss_grad(R1,R2,R3);
                diff2 = d * tau;
                error = norm(diff - diff2);
                tau = tau / 10;
                tau1 = tau1 / 10;
                tau2 = tau2 / 10;
                tau3 = tau3 / 10;
            end
        end
        function res = check_grad(this)
            tau1 = SO3.getRandomVector()';
            tau2 = SO3.getRandomVector()';
            tau3 = SO3.getRandomVector()';
            R1 = SO3(SO3.getRandomElement());
            R2 = SO3(SO3.getRandomElement());
            R3 = SO3(SO3.getRandomElement());
            tau = [tau1;tau2;tau3];
            tau = [tau1;tau2];
            for i = 1:9
                y = this.forward(R1,R2,R3);
                diff = this.forward(R1 * tau1, R2 * tau2, R3 * tau3) - y;
                d = this.grad(R1,R2,R3);
                diff2 = d * tau;
                error = norm(diff - diff2)
                tau = tau / 10;
                tau1 = tau1/10;
                tau2 = tau2/10;
            end
        end
        function res = check(this)
            tau1 = SO3.getRandomVector()';
            tau2 = SO3.getRandomVector()';
            R1 = SO3(SO3.getRandomElement());
            R2 = SO3(SO3.getRandomElement());
            tau = [tau1;tau2];
            for i = 1:9
            y = this.forward(R1,R2);
            diff = this.forward(R1 * tau1, R2 * tau2) - y;
            d = this.grad(R1,R2);
            diff2 = d * tau;
            error = norm(diff - diff2)
            tau = tau / 10;
            tau1 = tau1/10;
            tau2 = tau2/10;
            end
        end
        function res = check1(this)
            tau1 = [0;0;2];
            R1 = SO3(SO3.getRandomElement());
            tau = tau1;
            for i = 1:5
            y = R1.act(this.m_d1);
            new_R = (R1 * tau);
            diff = new_R.act(this.m_d1) - y;
            d = -Manifold.skew(R1.act(this.m_d1));
            diff2 = d * tau;
            error = norm(diff - diff2)
            tau = tau / 10;
            end
        end
        function res = check2(this)
            tau1 = SO3.getRandomElement();
            tau2 = SO3.getRandomElement();
            R1 = SO3(SO3.getRandomElement());
            R2 = SO3(SO3.getRandomElement());
            tau = [tau1;tau2];
            for i = 1:5
                y1 =  R1.act(this.m_d1);
                y2 =  R2.act(this.m_d2);
                y = y1 + y2;
                %y = y2;
                new_R1 = (R1 * tau1);
                new_R2 = (R2 * tau2);
                ny1 = new_R1.act(this.m_d1);
                ny2 = new_R2.act(this.m_d2);
                diff = ny1 + ny2 - y;
                %diff = ny2  - y;
                d1 = -Manifold.skew(R1.act(this.m_d1));
                d2 = -Manifold.skew(R2.act(this.m_d2));
                %d = [-Manifold.skew(R1.act(this.m_d1)),-Manifold.skew(R2.act(this.m_d2))];
                diff2 = [d1,d2] * [tau1;tau2];
                %diff2 = d2 * tau2;
                %diff2 = d * tau;
                error = norm(diff - diff2)
                tau1 = tau1/10;
                tau2 = tau2/10;
                tau = tau / 10;
            end
        end
        function d = grad(this, R1, R2, R3)
            grad_tau1 = - Manifold.skew(R1.act(this.m_d1)) ...
                        - R2.m_data * Manifold.skew(R1.act(this.m_d2)) ...
                        - R3.m_data * R2.m_data * Manifold.skew(R1.act(this.m_d3));
            grad_tau2 = - Manifold.skew(R2.act(R1.act(this.m_d2))) ...
                        - R3.m_data * Manifold.skew(R2.act(R1.act(this.m_d3)));
            grad_tau3 = - Manifold.skew(R3.act(R2.act(R1.act(this.m_d3))));
            d = [grad_tau1, grad_tau2, grad_tau3];
        end
        function d = loss_grad(this, R1, R2, R3)
            y = this.forward(R1, R2, R3);
            grad_manif = this.grad(R1, R2, R3);
            d = 2 * (y - this.m_target)' * grad_manif;
        end
    end
end