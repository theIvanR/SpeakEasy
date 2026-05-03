function plot_Ze(results)
%PLOT_ZE  Electrical impedance magnitude + phase.

    f  = results.f(:).';
    Ze = results.Ze;

    figure('Name','Electrical Impedance Z_e','Position',[400 400 1000 600]);

    subplot(2,1,1);
    semilogx(f, abs(Ze), 'LineWidth', 1.2);
    grid on;
    ylabel('|Z_e| (\Omega)');
    title('Electrical Impedance Magnitude');

    subplot(2,1,2);
    semilogx(f, unwrap(angle(Ze)) * 180/pi, 'LineWidth', 1.2);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Phase (degrees)');
    title('Electrical Impedance Phase');
end