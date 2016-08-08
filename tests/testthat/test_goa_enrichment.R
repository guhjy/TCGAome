context("Tests goa_enrichment.R functions")


## Test dataset for enrichment (source MCIA on TCGA's BRCAvsOV)
genes = c("ZNF638", "HNRNPU", "PPIAL4G", "RAPH1", "USP7", "SUMO1P3",
          "TMEM189.UBE2V1", "ZNF837", "LPCAT4", "ZFPL1", "STAT3", "XRCC1",
          "STMN1", "PGR", "RB1", "KDR", "YBX1", "YAP1", "FOXO3", "SYK", "RAB17",
          "TTC8", "SLC22A5", "C3orf18", "ANKRA2", "LBR", "B3GNT5", "ANP32E",
          "JOSD1", "ZNF695", "ESR1", "INPP4B", "PDK1", "TSC2", "AR", "HSPA1A",
          "CDH3", "SMAD4", "CASP7", "GMPS", "NDC80", "EZH2", "MELK", "CDC45",
          "CRY2", "KLHDC1", "MEIS3P1", "FBXL5", "EHD2", "CCNB1", "GSK3A",
          "DVL3", "NFKB1", "COL6A1", "CCND1", "BAK1")

entrez_ids = c("57763", "81611", "367", "84002", "578", "51161", "840", "891",
               "595", "8318", "1001", "1291", "1408", "1857", "30846", "2099",
               "2146", "26234", "2309", "8833", "2931", "3192", "3303", "8821",
               "9929", "3791", "122773", "3930", "254531", "9833", "10403",
               "4790", "5163", "5241", "644591", "64284", "65059", "5925",
               "6584", "4089", "6774", "3925", "100500808", "6850", "7249",
               "123016", "7874", "7515", "10413", "4904", "7542", "27332",
               "57116", "116412")

test_that("get_go_enrichment()", {
    ## Loads a GOA for a non supported universe
    results = TCGAome::get_go_enrichment(genes)
    testthat::expect_null(other_goa)
})
