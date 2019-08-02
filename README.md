# Deep learning models predict regulatory variants in pancreatic islets and refine type 2 diabetes association signals 

1. CNN training using scripts from Basset GitHub repository https://github.com/davek44/Basset

[Set up CNN training](CNN.create_input_and_training.sh) - creates the input files and runs the perl wrapper for training 100 nets with a specific param file
 
[Perl wrapper - CNN training](run_100_CNNs.pl) - emulates queuing system on GPU's, runs few CNN jobs at a time, waits and checks for available GPU resources to start the next CNN training
  
2. Evaluate CNN performance

[AUC per feature - per net](CNN.per_feature_AUC.Basset.sh)
  
[Validation set performance](validation_set.performance.sh)

[AUROC & AUPRC boxplots](AUROC_AUPRC.1000net_summary.R)
  
3. TF motif analysis

[Detect matches to known TF motifs](CNN.TF_motifs_analysis.sh)

[Motif annotation by filter size](motifs.R)

4. Regulatory predictions on HRC credible sets of variants

[Predict variant effects](CNN.HRC_credset_predictions.sh)

[Calculate CNN q-values](variant_CNN_q_values.per_feature.R)

[Append FGWAS and CNN scores to credible sets](credsets_add_FGWAS_CNN_Q.R)

[Convergence of fine-mapping PPAs and CNN predictions](convergence.R)

[Select signals better resolved by CNNs](CNN_fine_mapping.R)

5. Functional validation - ATAC-seq allelic imbalance

[Run ASEQ on islet ATAC-seq data](allelic_imbalance.ASEQ.sh)

[Allelic imbalance - binomial test](allelic_imbalance.binomial_test.R)

[Find islets heterozygous for PROX1 rs177](find_islets_het_for_PROX1_rs177.sh)

[Allelic imbalance at PROX1 signal](allelic_imbalance.PROX1.R)

[Enrichment of caQTLs among reg variants](enrichments_of_CNN_scores.incl_caQTLs.R)

6. Functional validation - enrichment of regulatory variants within islet chromHMM states

[ChromHMM states enrichment](chromHMM_enrichment.R)

7. Functional validation - PROX1 luciferase results plot

[Luciferase barplot](endoC_luciferase_plot.R)

8. Comparison with DeepSEA

[Compare to DeepSEA](compare_DeepSEA_and_my_CNN.081118.R)
