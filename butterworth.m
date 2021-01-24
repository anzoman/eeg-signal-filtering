% ------------------------------------------------------------------------
% 2. seminarska naloga pri predmetu KČR
% Avtor: Anže Luzar
% Tema: Izločanje očesnih artefaktov z uporabo postopka s pasovno 
%       prepustnim filtrom (Butterworth) in pragovno metodo
% ------------------------------------------------------------------------

function seminarska()
    %remove_eog_artifacts('S001R01.edf');
    
    % ostali primeri klicev funkcije:
    % remove_eog_artifacts('S001R01.edf', 'channel', 22);
    % remove_eog_artifacts('S001R01.edf', 'treshold', 75);
    % remove_eog_artifacts('S001R01.edf', 'channel', 22, 'treshold', 150, 'gap', 10);
    remove_eog_artifacts('S001R01.edf', 'channel', 22, 'treshold', 80, 'gap', 4, 'median', true);
end

% funkcija za odstranjevanje očesnih artefaktov
% sprejme en obvezen argument - pot do datoteke z EEG signalom in 4
% opcijske argumente - indeks kanala signala, prag filtriranja in velikost
% okna za pragovno metodo in še opcijski prikaz uprabe medianinega filtra
function remove_eog_artifacts(recordingPath, varargin)

  % privzete vrednosti za opcijske argumente
  options = struct('channel', 0, 'treshold', 150, 'gap', 10, 'median', false);

  % preberemo imena opcijskih argumentov
  optionNames = fieldnames(options);

  % preštejemo argumente
  nArgs = length(varargin);
  if round(nArgs/2)~=nArgs/2
    error('Opcijske argumente je potrebnopodati kot pare.')
  end

  for pair = reshape(varargin, 2, []) % par je {propName;propValue}
    inpName = lower(pair{1}); % case insensitive

    if any(strcmp(inpName, optionNames))
      % prepišemo privzete vrednosti opcijskih argumentov
      options.(inpName) = pair{2};
    else
      error('%s ni pravo ime parametra.', inpName)
    end
  end

   % iz datoteke s posetkom preberemo signal, frekvenco vzorčenja (kolikokrat 
   % na sekundo zajemamo signal) in časovno vrsto (časi zajema vzorcev)
   % če smo podali indeks kanala preberemo samo ta kanal signala, če ne pa
   % preberemo vseh prvih 64 kanalov signala (65. anotacijskega izpustimo)
  if options.('channel') ~= 0
    [eeg, freq, tm] = rdsamp(recordingPath, options.('channel'));
  else
    [eeg, freq, tm] = rdsamp(recordingPath, 1:64);
  end
  
  % normaliziramo frekvenco (delimo z 2)
  % uporabimo pasovno prepustni filter Butteworth (0.1 – 30 Hz) z n = 5
  [b1, a1] = butter(5, 0.1/(freq/2), 'high');
  [b2, a2] = butter(5, 30/(freq/2), 'low');
  % [b, a] = butter(5, [0.1 30]/(freq/2));

  % uporabimo dvosmerno shemo filtriranja
  eeg_butter_high = filtfilt(b1, a1, eeg);
  eeg_butter_low = filtfilt(b2, a2, eeg);
  eeg_butter_dvosmerno = filtfilt(b2, a2, eeg_butter_high);
  % eeg_butter_dvosmerno = filtfilt(b, a, eeg);

  % izvedemo pragovno metodo
  % nastavimo prag in N, ki pomeni koliko vzorcev skupaj odstranimo ko
  % najdemo prag (N/2 levo in N/2 desno)
  threshold = options.('treshold');
  N = options.('gap');
  eeg_butter_dvosmerno_tresholding = eeg_butter_dvosmerno;
  
  % poiščemo maksimume, ki presegajo prag in jih odstranimo (postavimo na
  % 0) ter v njihovi okolici odstranimo še N/2 elementov levo in N/2 desno
  idcs_max = find(eeg_butter_dvosmerno_tresholding > threshold);
  % dobljenim indeksom, ki jih je treba odstraniti dodamo še indekse, ki
  % so N/2 levo ali desno od določenega indeksa in dobimo novo matriko
  idcs_max_n = unique(cell2mat(arrayfun(@(x) x-(N/2):x+(N/2), idcs_max, 'UniformOutput', false)));
  idcs_max_n = idcs_max_n(idcs_max_n <= max(idcs_max));
  eeg_butter_dvosmerno_tresholding(idcs_max_n) = 0;
  
  % poiščemo minimume, ki so manjši od -prag in jih odstranimo ter v 
  % njihovi okolici odstranimo še N/2 elementov levo in N/2 desno
  icds_min = find(eeg_butter_dvosmerno_tresholding < -threshold);
  icds_min_n = unique(cell2mat(arrayfun(@(x) x-(N/2):x+(N/2), icds_min, 'UniformOutput', false)));
  icds_min_n = icds_min_n(icds_min_n >= min(icds_min));
  eeg_butter_dvosmerno_tresholding(icds_min_n) = 0;
  
  % narišemo grafe
  figure(1);
  if options.('median')
    tiledlayout(4, 1);
  else
    tiledlayout(3, 1);
  end
  % prvi podgraf
  ax1 = nexttile;
  plot(tm, eeg);
  title('Originalni signal');
  xlabel('Časovni indeks');
  ylabel('Signal');
  % drugi podgraf
  ax2 = nexttile;
  plot(tm, eeg_butter_dvosmerno);
  title('Signal po dvosmerni obdelavi z Butteworth filtrom');
  xlabel('Časovni indeks'); 
  ylabel('Signal');
  % tretji podgraf
  ax3 = nexttile;
  plot(tm, eeg_butter_dvosmerno_tresholding);
  title('Signal po obdelavi s pragovno metodo');
  xlabel('Časovni indeks');
  ylabel('Signal');
  
  if options.('median')
    % četrti podgraf (narišemo signal po obdelavi z 1D medianinim filtrom)
    ax4 = nexttile;
    plot(tm, medfilt1(eeg_butter_dvosmerno, 100));
    title('Signal po obdelavi z 1D medianinim filtrom');
    xlabel('Časovni indeks');
    ylabel('Signal');
    
    % uskladimo osi podgrafov
    linkaxes([ax1 ax2 ax3 ax4], 'xy');
  else  
    % uskladimo osi podgrafov
    linkaxes([ax1 ax2 ax3], 'xy');
  end
  
  % če smo izbrali vse kanale narišemo vsakega posebej pred in po operaciji
  % dvosmernega filtriranja z Butterworth filtrom in pragovno metodo
  if options.('channel') == 0
    figure(2);
    tl_1 = tiledlayout(8, 8);
    axes = [];
    for i=1:size(eeg, 2)-1
      axes = [axes, nexttile(tl_1)];
      plot(tm, eeg(:, i));
      hold on;
      plot(tm, eeg_butter_dvosmerno_tresholding(:, i));
      hold off;
    end
    title(tl_1, 'Vseh 64 originalnih kanalov signala in obdelani dvosmerno s filtrom Butterwoth in s pragovno metodo');
    linkaxes(axes, 'xy');
  end

end
    