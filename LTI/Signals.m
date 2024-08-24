classdef Signals < Scene
  properties (Access = private)
    m_signals Signal
    m_func
    m_temp
  end

  methods
    function this = Signals()
      this@Scene("signals.gif", "function");
      this.m_func = @(t) 0;
    end

    function this = plus(this, signal)
      this.m_signals = [this.m_signals signal];
      func = @(t) 0;
      for i = 1:size(this.m_signals,2)
          temp = func;
          func = @(t) temp(t) + this.m_signals(i).m_func(t);
      end
      %debug
      func(1);
      %enddebug
      this.m_func = func;
    end

    function Show(this)
        n = 1000;
        T = 10;
        t = linspace(0,T,n);
        y = this.m_func(t);
        plot(real(y),imag(y));
    end

    function ShowDynamically(this)
        n = 200;
        T = 20;
        t = linspace(0,T,n);
        y = this.m_func(t);
        this.setBB([real(y); imag(y)]);
        for i = 1:n
            h = [];
            y_i = [0,0];
            for j = 1:size(this.m_signals,2)
                s_i = this.m_signals(j).m_func(t(i));         
                s_i = [real(s_i), imag(s_i)];
                h = [h,quiver(y_i(1), y_i(2), s_i(1), s_i(2), 'color', this.m_red, 'AutoScale','off')];
                y_i = y_i + s_i;
            end
            y_j = y(1:i);
            h = [plot(real(y_j),imag(y_j),'color',this.m_blue, 'linewidth',1.0),h];
            this.get();
            delete(h);
        end
        this.save();
    end
  end
end