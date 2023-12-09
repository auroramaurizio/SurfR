\name{enrichedList}
\alias{enrichedList}
\title{enrichedList}
\value{
  A list of lists.
}
\description{
  Input list for Enrichment_barplot function. enrichedList is the output of
  Enrichment function applied to ind_deg object when enrich.databases is equal
  to GO_Cellular_Component_2021, default parameters.
}
\usage{
  data(enrichedList)
}
\format{
  list
}
\details{
  enrichedList$fdr_up$GO_Cellular_Component_2021 contains upregulated gene
  enrichments, enrichedList$fdr_down$GO_Cellular_Component_2021 contains
  downregulated gene enrichments.
}

