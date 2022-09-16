# 2022.09.16 Skype meeting (Dave and Veronica):

(A) MB HIJING ("A and B") --> provides background smearing for FastJet clustered jets. DONE
(B) PYTHIA embedded in HIJING --> This gives us actual response Truth::fulljets to Detector::calojets. DONE
(C) PYTHIA not embedded --> Provide Truth to Detector w/o background. NEEDED

Will do: 
Take (C: TH2D) Response matrix per centrality bin, smear with (A: TH1D)
- Compare to (B) for how good it is (and compare to B's closure)
- Do closure and see how good closure it
- [ ] Do (C) : Veronica will run (B) with input files that don't have HIJING background -- are just PYTHIA jets
