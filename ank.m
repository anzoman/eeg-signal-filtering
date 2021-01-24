% ------------------------------------------------------------------------
% Izpitna seminarska naloga pri predmetu KČR
% Avtor: Anže Luzar
% Tema: Izločanje očesnih artefaktov z uporabo postopka 
%       analize neodvisnih komponent (ANK)
% ------------------------------------------------------------------------

function izpit()
    remove_eog_artifacts('S010R01.edf', 'remove', [2, 3, 10, 19, 23, 35, 36]);
    
    % ostali primeri klicev glavne funkcije:
    % remove_eog_artifacts('S010R01.edf');
    % remove_eog_artifacts('S001R01.edf', 'channel', 24);
    % remove_eog_artifacts('S001R01.edf', 'remove', [1, 2, 3, 4, 5]);
    % remove_eog_artifacts('S001R01.edf', 'topography', true);
    % remove_eog_artifacts('S001R01.edf', 'channel', 24, 'remove', [1, 2, 3], topography', true);
end

% funkcija za odstranjevanje očesnih artefaktov
% sprejme en obvezen argument - pot do datoteke z EEG signalom in 3
% poimenovane opcijske argumente - indeks kanala signala, matriko, ki 
% vsebuje naštete ANK komponente za izločanje in argument za topografijo,
% ki izriše dodatne topografske slike posameznih komponent signala
function remove_eog_artifacts(recordingPath, varargin)

  % privzete vrednosti za opcijske argumente
  options = struct('channel', 0, 'remove', [], 'topography', false);

  % preberemo imena opcijskih argumentov
  optionNames = fieldnames(options);

  % preštejemo argumente
  nArgs = length(varargin);
  if round(nArgs/2)~=nArgs/2
    error('Opcijske argumente je potrebno podati kot pare.')
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

   % iz datoteke s posnetkom preberemo signal, frekvenco vzorčenja (kolikokrat 
   % na sekundo zajemamo signal) in časovno vrsto (časi zajema vzorcev)
   % preberemo vseh prvih 64 kanalov signala (65. anotacijskega izpustimo)
   [eeg, freq, tm] = rdsamp(recordingPath, 1:64);
  
  % narišemo 64 osnovnih kanalov prebranega EEG signala
  figure('Name', 'Osnovni EEG signali');
  tl_1 = tiledlayout(8, 8);
  axes = [];
  for i=1:size(eeg, 2) - 1
    axes = [axes, nexttile(tl_1)];
    plot(tm, eeg(:, i));
  end
  title(tl_1, 'Vseh 64 originalnih kanalov signala');
  linkaxes(axes, 'xy');
  
  % obrnemo matriko s signalom, saj fastica pričakuje signale po vrsticah
  sig = eeg';
  
  % uporabimo fastica skripto da zeženemo metodo ANK nad vhodnim EEG signalom 
  % dobimo neodvisne komponente icasig, mešalno matriko A in ločevalno matriko W
  [icasig, A, W] = fastica(sig);
  
  % narišemo signale EEG v prostoru komponent dobljene z metodo ANK
  figure('Name', 'Signali EEG v prostoru komponent');
  tl_2 = tiledlayout(8, 8);
  axes = [];
  for i=1:size(icasig, 1) - 1
    axes = [axes, nexttile(tl_2)];
    plot(tm, icasig(i, :));
  end
  title(tl_2, 'Signali EEG v prostoru komponent');
  linkaxes(axes, 'xy');
  
  % izračunamo inverz ločevalne matrike W
  % dobimo novo mešalno matriko W1
  W1 = inv(W);

  % če uporabnik želi izrišemo topografske slike signalov komponent
  % seznam vseh 64 kanalov 
  ch_list = {'Fc5', 'Fc3', 'Fc1', 'Fcz', 'Fc2', 'Fc4', 'Fc6', ...
           'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', ...
           'Cp5', 'Cp3', 'Cp1', 'Cpz', 'Cp2', 'Cp4', 'Cp6', ...
           'Fp1', 'Fpz', 'Fp2', ...
           'Af7', 'Af3', 'Afz', 'Af4', 'Af8', ...
           'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8', ...
           'Ft7', 'Ft8', ...
           'T9', 'T7', 'T8', 'T10', ...
           'Tp7', 'Tp8', ...
           'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', ...
           'Po7', 'Po3', 'Poz', 'Po4', 'Po8', ...
           'O1', 'Oz', 'O2', ...
           'Iz'};

  if options.('topography')
    figure('Name', 'Topološka slika signalov');
    j = 0;
    for i = 1:64
      j = j + 1;
      subplot(4, 4, j);
      plot_topography(ch_list, W1(:, i));
      title(sprintf('Komponenta %d', i))
    
      if mod(i, 16) == 0
        figure('Name', 'Topološka slika signalov');
        j = 0;
      end
    end
  end

  % uporabnik izbere katere komponente bo izločil, če jih ni podal že prej
  if isempty(options.('remove'))
    prompt = 'V matriki naštej katere komponente želiš izločiti: ';
    izloci = input(prompt)
    if ~ismatrix(izloci)
      exit
    end
  else
    izloci = options.('remove');
  end
  
  % iz signala izbrišemo neodvisne komponente, ki vsebujejo šum 
  % odstranimo ustrezne stolpce iz matrike W1 in dobimo mešalno matriko Wap
  Wap = W1;
  Wap(:, izloci) = [];
  Yap = icasig
  Yap(izloci, :) = [];

  % izračunamo in narišemo nov korigiran signal brez neodvisnih komponent
  Xap = Wap * Yap;
  figure('Name', 'Korigirani EEG signali');
  tl_3 = tiledlayout(8, 8);
  axes = [];
  for i=1:size(Xap, 1) - 1
    axes = [axes, nexttile(tl_3)];
    plot(tm, Xap(i, :));
  end
  title(tl_3, 'Korigirani EEG signali');
  linkaxes(axes, 'xy'); 
  
  % narišemo vse kanale vsakega posebej pred in po operaciji
  % postopka analize neodvisnih komponent (ANK)
  figure('Name', 'Vseh 64 originalnih kanalov signala in obdelani signali z metodo ANK');
  tl_4 = tiledlayout(8, 8);
  axes = [];
  for i=1:size(eeg, 2)-1
    axes = [axes, nexttile(tl_4)];
    plot(tm, eeg(:, i));
    hold on;
    plot(tm, Xap(i, :));
    hold off;
  end
  title(tl_4, 'Vseh 64 originalnih kanalov signala in obdelani signali z metodo ANK');
  linkaxes(axes, 'xy');
  
  % narišemo preostale grafe
  figure('Name', 'Rezultat filtriranja');
  tiledlayout(2, 1);
  % prvi podgraf
  ax1 = nexttile;
  if options.('channel') ~= 0
      plot(tm, eeg(:, options.('channel')));
  else
      plot(tm, eeg);
  end
  title('Originalni vhodni EEG signal');
  xlabel('Časovni indeks');
  ylabel('Signal');
  % drugi podgraf
  ax2 = nexttile;
  if options.('channel') ~= 0
    plot(tm, Xap(options.('channel'), :));
  else
    plot(tm, Xap);
  end
  title('EEG signal po odstranjevanju očesnih artefaktov z metodo ANK');
  xlabel('Časovni indeks'); 
  ylabel('Signal');
  % uskladimo osi podgrafov
  linkaxes([ax1 ax2], 'xy');
  
end
 