function [nCPU, totCPU, nBlockPerCPU, totSLAVES] = distributeJobs(Parallel, fBlock, nBlock)
% PARALLEL CONTEXT
% In parallel context this function is used to determine the total number of available CPUs,
% and the number of threads to run on each CPU.
%
% INPUTS
%  o Parallel [struct vector]   copy of options_.parallel
%  o fBlock [int]               index number of the first job (e.g. MC iteration or MH block)
%                               (between 1 and nBlock)
%  o nBlock [int]               index number of the last job.
%
% OUTPUT
%  o nBlockPerCPU [int vector]  for each CPU used, indicates the number of
%                               threads run on that CPU
%  o totCPU [int]               total number of CPU used (can be lower than
%                               the number of CPU declared in "Parallel", if
%                               the number of required threads is lower!)
%  o nCPU                       the number of CPU in user format.
%  o totSLAVES                  dovrebbe rappresentare il numero dei nodi
%                               di calcolo elencati in Parallel ed
%                               effettivamente coinvolti nella computazione
%                               attuale è compreso tra 1 e
%                               length(Parallel).


% Copyright (C) 2010 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

totCPU=0;
for j=1:length(Parallel),
    nCPU(j)=length(Parallel(j).CPUnbr);
    totCPU=totCPU+nCPU(j);
end

nCPUoriginal=nCPU;

nCPU=cumsum(nCPU);
offset0 = fBlock-1;
if (nBlock-offset0)>totCPU,
    diff = mod((nBlock-offset0),totCPU);
    nBlockPerCPU(1:diff) = ceil((nBlock-offset0)/totCPU);
    nBlockPerCPU(diff+1:totCPU) = floor((nBlock-offset0)/totCPU);
    totSLAVES=length(Parallel);
else
    nBlockPerCPU(1:nBlock-offset0)=1;
    totCPU = nBlock-offset0;
    totSLAVES = min(find(cumsum(nCPU)>=totCPU));
end

% Supponiamo che ereditiamo un vettore normalizzato
% della lunghezza di Parallel con tutti i valori > 0.
% Per avere un valore 0 basta non elencarlo sopra nei nodi coinvolti nel
% calcolo o non mettere il  nodo nel file di configurazione.
% Supponiamo inoltre che tutti i controlli per avere questa consistenza
% siano fatte dal compilatore o nella Analyze ...

% La notra filosofia fino ad ora è:
% 1.    Considera la mole di lavoro che devi fare,
% 2.    Valuta le risorse che hai,
% 3.    Parti dalla risorsa numero 1 e distribuisci il lavoro tra le
%       diverse risorse in modo bilanciato.

% Con questa soluzione tutte le cpu, partendo dal primo nodo elencato nel
% file di configurazione hanno lo stesso carico. L'unica eccezione può
% essere l'ultima cpu quando le configuazioni possibili sono 'dispari' e
% viene a trovarsi con un job in meno.

% Io per modificare il meno possibile farei semplicemente così:

% Faccio i punti 1, 2, 3 come prima, e in questo punto ho:


% nCPU:         numero delle cpu dichiarate nel file di configurazione,

% totCPU:       numero delle cpu che verranno effettivamente
%               coinvolte nello step parallelo considerato. Ogni volta viene ricalcolato
%               tutto.
%
% nBlockPerCPU  è il numero di threads che deve eseguire una cpu/core (non il
%               nodo che può avere molte cpu e molti core).

% totSLAVES     numero delle macchine (non cpu) effettivamente coinvolte
%               nell'attuale frazione di calcolo parallelo.


% Esempi:

% 1.    mh_nblocks=2,mh_replic=1005
%       con 2 nodi il primo con 2 cpu il secondo con 1 cpu
%
%       Durante la computazione possiamo avere situazioni come:
%
%       fBlock =1
%       nBlock =2
%       nCPU = 2     3
%       totCPU = 2
%       nBlockPerCPU = 1     1
%       totSlaves = 1

%       Che significa:
%       devo fare 2 jobs ho due macchine la prima con 2 cpu
%       la terza con 3-2=1 cpu.
%
%       Quindi per farli uso solo la prima e gli assegno un job a cpu.
%       Il secondo nodo è inattivo.
%
%
%
% 2. Allo stesso modo se ho:
%    fBlock =1
%    nBlock = 17
%    nCPU = 2     3
%    totCPU = 3
%    nBlockPerCPU = 6     6     5
%    totSlaves = 2
%
% Significa:
%       devo fare 17 jobs ho due macchine la prima con 2 cpu
%       la terza con 3-2=1 cpu.
%
%       Per farli le uso tutte e tre e assegno 6 job alla prima 6 alla seconda
%       e 5 all'ultima.
%       Tutti e due i nodi sono attivi.


% Quindi per evitare di cambiare il codice posso semplicemente fare cosi:
%
% Abbiamo
% lp= length(Parallel);
% CPUWeight[c1 c2 ... clp];
% Con c1+c2+...+clp=1;
%
% lc=length(nCpu);

% if (Tutti i ci sono uguali) | (L'utente non definisce CPUWeight)
%
%           NON FARE NIENTE perchè:
%           - Quello di adesso va bene.
%             oppure perchè
%           - Mi viene chiesto di non fare niente.
%
% else
%     Considera tutti i nodi,
%     Per tutti i ci in CPUWeight, fai:
%     Considera la frazione ci del numero totale dei jobs (=
%     nBlock-fBlock+1) e assegnali al nodo ni. Se il nodo ni
%     ha più di una cpu, distribuiscili in modo uniforme tra le cpu.
% end

%   Possibile Implementazione

global options_

% Copio in locale e normalizzo ...
CPUWeight=options_.CPUWeight.*nCPUoriginal,
CPUWeight=CPUWeight/sum(CPUWeight)


lCw=length(CPUWeight);

EqFlag=1;

for i=1:(lCw-1)
    if CPUWeight(i)~=CPUWeight(i+1)
        EqFlag=0;
        SonoQui='Diverso'
        break;
    end
end

% L'utente non ha inserito il vettore di pesi, oppure i pesi sono tutti
% uguali.

if (EqFlag==1) | (lCw==0)
    SonoQui='Uguale'
    return;
    
else
    
    % Numero dei Nodi nel cluster ...
    lnC=length(nCPUoriginal);
    
    % Numero totale dei Jobs ...
    NumbersOfJobs=sum(nBlockPerCPU);
    
    
    SumOfJobs=0;
    JobsForNode=zeros(1,lnC);
    
  %    keyboard
    
    % Ridistribusco i jobs tra i nodi in base ai pesi dell'utenti.
    
    for i=1:lnC
              
        JobsForNode(i)=CPUWeight(i)*NumbersOfJobs;        
        % Ci sono diverse soluzioni possibili: round sembra la
        % migliore.
        
        JobsForNode(i)=ceil(JobsForNode(i));
%         JobsForNode(i)=round(JobsForNode(i));
        
      %  SumOfJobs=sum(JobsForNode(1:i));
    end     
        
        % Tolgo gli eventuali 'eccessi' o 'mancanze' derivanti dall modo usato sopra
        % per 'arrotondare'.
        % Ci sono diverse strategie per fare questo ...
        
          SumOfJobs=sum(JobsForNode);
          
        if SumOfJobs~=NumbersOfJobs
            
            
            % Ho assegnato più jobs di quelli veri ..
            
            if SumOfJobs>NumbersOfJobs
                % Li tolgo al meno veloce,
                % posso anche toglierli a chi ne ha di più, ... da decidere
                % ...
                
                [NonServe VerySlow]= min(CPUWeight);
               
                while SumOfJobs>NumbersOfJobs
                    JobsForNode(VerySlow)=JobsForNode(VerySlow)-1;
                    SumOfJobs=SumOfJobs-1;
                end
                
            end
            
            if SumOfJobs<NumbersOfJobs
                
                % Li metto al più veloce,
                % posso anche toglierli a chi ne ha di più, ... da decidere
                % ...
                
                [NonServe VeryFast]= min(CPUWeight);
                
                while SumOfJobs<NumbersOfJobs
                    JobsForNode(VeryFast)=JobsForNode(VeryFast)+1;
                    SumOfJobs=SumOfJobs+1;
                end
                
            end
        end
        
   
    
    
    % Adesso ridistribusco i jobs assegnati ad ogni nodo tra le
    % cpu/core disponibili in quel nodo! Poi si può eventualmente
    % accorpare con il codice sopra.
    
    JobsForCpu=zeros(1,nCPU(lnC));
    JobAssignedCpu=0;
    
    RelativePosition=1;
    
    
    
    for i=1:lnC
        
        
        % Diverse possibilità ...
        
        % JobAssignedCpu=ceil(JobsForNode(i)/nCPUoriginal(i));
%         JobAssignedCpu=round(JobsForNode(i)/nCPUoriginal(i));
        JobAssignedCpu=floor(JobsForNode(i)/nCPUoriginal(i));
        
        ChekOverFlow=0;
        
        for j=RelativePosition:nCPU(i)
            JobsForCpu(j)=JobAssignedCpu;
            ChekOverFlow=ChekOverFlow+JobAssignedCpu;
            
            if ChekOverFlow>=JobsForNode(i)
                break;
            end
            
        end
        
        % Tolgo gli eventuali 'eccessi'- 'mancanze' derivanti dall modo usato sopra
        % per 'arrotondare'. Anche qui come sopra è da decidere la
        % strategia migliore ...
        
        if ChekOverFlow ~=(JobsForNode(i))
            
            if ChekOverFlow >(JobsForNode(i))
                while ChekOverFlow>JobsForNode(i)
                    JobsForCpu(nCPU(i))=JobsForCpu(nCPU(i))-1;
                    ChekOverFlow=ChekOverFlow-1;
                end
            end
                 
            if ChekOverFlow <(JobsForNode(i))
                while ChekOverFlow<JobsForNode(i)
                    JobsForCpu(nCPU(i))=JobsForCpu(nCPU(i))+1;
                    ChekOverFlow=ChekOverFlow+1;
                end
            end
        end
        
        RelativePosition=nCPU(i)+1;
        
    end
    
    % Only for testing ...
    
    SonoQui='Maggiore'
    
    nCPUoriginal
    nCPU
    JobsForNode
    
    
    display('New')
    
    %nBlockPerCPUWweigh= JobsForCpu(JobsForCpu~=0)
    JobsForCpu
    
    display('Old')
    
    nBlockPerCPU
    
    display('Check')
    
    Check=sum(JobsForCpu)-sum(nBlockPerCPU)
    
    %pause
end


