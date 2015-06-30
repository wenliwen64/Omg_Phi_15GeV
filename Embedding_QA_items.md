## Embedding QA items

The embedding process is: embed the known, Monte-Carlo tracks of exact particle type to detector geometry and run through real track/event reconstruction procedure to get embedding data files.
###$$\Omega^\pm$$ for Run12 AuAu 14.5GeV
###Pt Bin
---
Bin no. | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 
------- |---|---|---|---|---|---|---|
Pt(GeV/c)|(0.7, 1.2]|(1.2, 1.6]|(1.6, 2.0]|(2.0, 2.4]|(2.4, 2.8]|(2.8, 3.6]|(3.6, $$\inf$$)

###Cuts
---
- Event-wise:
   - [X] VertexZ;
   - [X] VertexR;
   - Bad runs? (How embedding works? Can it simulate detector performance on daily/run basis?)
 
- Track-wise:
   - Hits number;
   - Pt; 

- Lambda-wise:
   - proton_dca;
   - pion_dca;
   - dca\_proton\_to_pion;
   - lambda_dca;
   - lambda\_decay_length;
   - |lambda_mass - PDG|;

- Omega-wise:
   - bach_dca;
   - dca\_bach_to\_v0;
   - omega_dca;
   - omega_decay\_length;
   - omega_decay\_length(< lambda\_decay\_length);
   - omega_rapidity;
   - omega_mass within pdgmass$$\pm$$0.008GeV(when calculating spectra);

- Colinear cuts:
   - $$(\vec{r}_\Omega - \vec{r}_\Lambda)\cdot \vec{p}_\Lambda$$ > 0
   - $$(\vec{r}_\Omega - \vec{r}_{pv})\cdot \vec{p}_\Omega$$ > 0
   - $$\frac{|(\vec{r}_\Omega - \vec{r}_{PV})\times \vec{p}_\Omega|}{|\vec{r}_\Omega - \vec{r}_{PV}|\cdot|\vec{p}_\Omega|}$$  < 0.15
