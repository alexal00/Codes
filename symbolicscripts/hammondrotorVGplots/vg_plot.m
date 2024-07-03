function out=vg_plot(p,U)
    sz = size(p,1);
    lambda = p.*U;
    ax2 = subplot(2,1,2);
    % Use g as damping factor
    % plot(U*60/(2*pi), 2*real(lambda(1,:))./(abs(imag(lambda(1,:)))+1*(abs(imag(lambda(1,:)))< 0.001)), '*-');
    % Use zeta as damping factor
    plot(U*60/(2*pi), real(lambda(1,:)), '*-');
    xlabel(ax2, '$\Omega$ [rpm]');
    %ylabel('Damping ratio, %  ');
    ylabel(ax2, '$\Re (\lambda)$ [rad/s]');
    hold on;
    for i = 2:sz
        % plot(U*60/(2*pi), 2*real(lambda(i,:))./(abs(imag(lambda(i,:)))+1*(abs(imag(lambda(i,:)))< 0.001)), '*-');grid on;
        plot(U*60/(2*pi), real(lambda(i,:))./(sqrt(imag(lambda(1,:)).^2+real(lambda(i,:)).^2)), '*-');
    end
    grid on;
    
    ax1 = subplot(2,1,1);
    plot(U*60/(2*pi), abs(imag(lambda(1,:)))/2/pi, '*-');
    xlabel(ax1, '$\Omega$ [rpm]');
    ylabel(ax1, '$\Im (\lambda)$ [rad/s]');
    hold on;
    for i = 2:sz
        plot(U*60/(2*pi), abs(imag(lambda(i,:)))/2/pi, '*-');
    end
    grid on;
end