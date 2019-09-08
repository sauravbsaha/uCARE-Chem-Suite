#Loading the library
#Software Name: uCARE Chem Suit Version 0.2
#Authors: Saurav Bhaskar Saha, Pramod Wasudeo Ramteke
#Details: The script takes sdf file as input and produces output through
#		  several tabs viz. Drug Chemical Properties,

library("shiny")
library("ChemmineR")
library("ape")
library("uCAREChemSuiteCLI")

rm(list = ls())
file_plot<- "www/query_db_cluster.pdf"
file_output_text<- "www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt"
file_out_pdf<- "www/uCARE.Chem.Suit.Version.0.2.Chemical.Structure.Output.pdf"

if (file.exists(file_plot)) file.remove(file_plot)
if (file.exists(file_output_text)) file.remove(file_output_text)
if (file.exists(file_out_pdf)) file.remove(file_out_pdf)

############################# UI File #########################
# UI Declaration
ui = fluidPage(
  headerPanel("uCARE Chem Suite V.02"),
  #Title of application
  titlePanel(h2("")
             #h2("uCARE Chem Suit V.02 ", align="center"), tags$img(src='images/logo1.png', height=125, width=285)
  ),

  #Side Bar Layout
  sidebarLayout(

    #Side Panel asking for scientist and Project details
    sidebarPanel(

      fileInput('file', "Choose drug SDF file to start project",
                accept='.sdf', placeholder = "No file selected"),
      h4(strong("Algorithm parameters")),
      radioButtons("organism","Pathogen name",list("Escherichia coli","Pseudomonas aeruginosa")),
      radioButtons("application","Application Type",list("Visualization", "Resistome Prediction")),
      radioButtons("model","Model Type",list("Deterministic Model", "Stochastic Model")),
      radioButtons("neighbor","Nearest Neighbor",list("1", "3")),
      radioButtons("threshold","Threshold structural similarity score",list("0.25", "0.30","0.35","0.40")),
      h5(strong("Project details")),
      textInput("projID","Enter the Project ID (Optional)", value = "NA"),
      textInput("name","Enter your name (Optional)", value = "NA"),
      textInput("email","Enter your Institutional email id (Optional)", value = 'NA'),
      textAreaInput("com","Comments/Notes (Optional)", value = "NA"),
      ("Note: Project details are optional headers used as description to exported output file. Moreover, for every new session, it deletes the history. "),br(),


      "Date:", date(),br()
    ),

    #Main Panel to print
    mainPanel(uiOutput('home') )
  )
)

############################# Server File #########################
server = function(input, output) {

  #Extracting Molecular Formula of drug SDF file
  output$mf = renderTable({
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    data.frame(MF=MF(read.SDFset(inFile$datapath)[[1]]), MW=MW(read.SDFset(inFile$datapath)[[1]]))
  },
  include.rownames=TRUE)

  ################################################################
  #Extraction of model Type
  output$mymodel <- renderText(input$model)

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #Extracting Atomblock of Drug SDF File
  output$mytable = renderTable({
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    atomblock(read.SDFset(inFile$datapath)[[1]])
  },
  include.rownames=TRUE, include.colnames=TRUE)

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #Extracting Bond block of Drug SDF File
  output$bondblock = renderTable({
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    bondblock(read.SDFset(inFile$datapath)[[1]])
  },
  include.rownames=TRUE, include.colnames=TRUE)

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #Extracting Ring attribute of Drug SDF File
  output$ring = renderTable({
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    rings(read.SDFset(inFile$datapath)[[1]],type="count",upper=1000, arom=TRUE)
  },
  include.rownames=TRUE, include.colnames=FALSE)

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #Extracting Bond matrix of Drug SDF File
  output$bmatrix = renderTable({
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    conMA(read.SDFset(inFile$datapath)[[1]], exclude=c("H"))
  },
  include.rownames=TRUE,include.colnames=TRUE)

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #Extracting Ring structure attribute of Drug SDF File
  output$rstructure = renderTable({
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    rc <- rings(read.SDFset(inFile$datapath)[[1]], upper=Inf, type="all", arom=TRUE, inner=FALSE)

    rc_df<- c(rc$RINGS[seq(from= 1, to = length(rc$RINGS))])
    rc_df <- as.character(rc_df)
    rc_df[rc_df==""] <- "NA"
    rc_df <- as.factor(rc_df)
    rc_df
  },
  include.rownames=TRUE, include.colnames=FALSE)

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #Extracting Number of Ring in Drug SDF File
  output$rnum = renderTable({
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    rc <- rings(read.SDFset(inFile$datapath)[[1]], upper=Inf, type="all", arom=TRUE, inner=FALSE)
    rc[2]
  },
  include.rownames=TRUE)

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #Extracting chemical structure of drug SDF file
  output$img = renderPlot({
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    plot(read.SDFset(inFile$datapath)[[1]], print=FALSE, atomnum=TRUE, no_print_atoms="H")
  })

  #########################################################################
  ###$$$$$$$$$$$$$$$$$$$$$$$$$ Rule based analysis $$$$$$$$$$$$$$$$$$$$$###
  #########################################################################
  #Nearest neighbor-Deterministic algorithm
  #Extracting Drug neighbors
  output$dclass = renderTable({
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    sdf_input<-read.SDFset(inFile$datapath)[[1]]
    write.SDF(sdf_input,"input.sdf")
    dc<-drug.class.deterministic("input.sdf")

    organism <- input$organism
    if(organism == "Escherichia coli")
    {
    db<- read.csv2("db/Antibiotic_data_set_Ecoli.csv", header = TRUE,sep=",")
    resistome_deterministic = data.frame(subset(db, Drug_class == dc))

    if(grepl("The drug is not found in the list!", dc) == TRUE)
    {print(paste("Deterministic model could not predict the class of input compound! Please use stochastic model."))}
      else
        {
          if(nrow(resistome_deterministic) == 0)
          {print(paste("No resistant drug for the predicted class could be found in the database for the specified genome."))}
          else{
            drug.neighbor<-data.frame(resistome_deterministic[,c(1,6)])
            subset(drug.neighbor, !duplicated(Drug_Name))
          }
        }
    }else{
      db<- read.csv2("db/Antibiotic_data_set_Paeruginosa.csv",header =TRUE, sep=",")
      resistome_deterministic = data.frame(subset(db, Drug_class == dc))

      if(grepl("The drug is not found in the list!", dc) == TRUE)
      {print(paste("Deterministic model could not predict the class of input compound! Please use stochastic model."))}
      else
        {
          if(nrow(resistome_deterministic) == 0){
            print(paste("No resistant drug for the predicted class could be found in the database for the specified genome"))
          }
          else{
            drug.neighbor<-data.frame(resistome_deterministic[,c(1,6)])
            subset(drug.neighbor, !duplicated(Drug_Name))
          }
        }
    }



  },
  include.rownames=FALSE, include.colnames=FALSE)

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$ Prediction Rules Ends here $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


  ########################################################################
  ################################ Drug Class classification #############
  ########################################################################

  # Extracting Drug class
  output$cclass = renderTable({
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    sdf_input<-read.SDFset(inFile$datapath)[[1]]
    write.SDF(sdf_input,"input.sdf")
    dc<-drug.class.deterministic("input.sdf")
   if(grepl("The drug is not found in the list!", dc) == TRUE)
   {print(paste("Deterministic model could not predict the class of input compound! Please use stochastic model."))}
    else{print(paste("The drug is likely to belong ",dc,"class."))}
  },
  include.rownames=FALSE, include.colnames=FALSE)


  ############################ Drug classification ends here##############

  #########################################################################
  ###$$$$$$$$$$$$$$$$$$$$$$$$$ Rule based Gene classification$$$$$$$$$$$$$$$$$$$$$###
  #########################################################################

  output$gclass = renderDataTable({

    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    sdf_input<-read.SDFset(inFile$datapath)[[1]]
    write.SDF(sdf_input,"input.sdf")
    organism <- input$organism
    if(organism == "Escherichia coli")
    {
      resistome_deterministic<- drug.resistome.deterministic("input.sdf", 1)
    }
    else
    {
      resistome_deterministic<- drug.resistome.deterministic("input.sdf", 2)
    }





  },escape = FALSE)

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$ Drug-Gene classification Ends here $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

  #######################################################################
  ######################## Export Tab ##################################
  #######################################################################

  checkbox_Drug_Chemical<- reactive({
    checkboxInput("Drug_Chemical_Properties", "Drug Chemical Properties", value = FALSE, width = NULL)

  })

  checkbox_Drug_Atomic<- reactive({
    checkboxInput("Drug_Atomic_Properties", "Drug Atomic Properties", value = FALSE, width = NULL)
  })

  checkbox_Bond_Attributes<- reactive({
    checkboxInput("Bond_Attributes", "Bond Attribute", value = FALSE, width = NULL)
  })

  checkbox_Drug_Classification<- reactive({
    checkboxInput("Drug_Classification", "Drug Classification", value = FALSE, width = NULL)
  })

  checkbox_Nearest_drug<- reactive({
    checkboxInput("Nearest_drugs", "Nearest drug/s", value = FALSE, width = NULL)
  })

  checkbox_Resistance_Gene_List<- reactive({
    checkboxInput("Resistance_Gene_List", "Resistance Gene List", value = FALSE, width = NULL)
  })



  #######################################################################################################################
  #######################################################################################################################

  ############################################### Export Functions ######################################################

  #######################################################################################################################
  #######################################################################################################################


  ######################################## View  and Download Files ###############################################
  export_submit<- reactive({
    actionButton("buttonuno","View Output")
  })

  export_download<- reactive({
    actionButton("buttondos","Create File")

  })



  values <- reactiveValues(uno =0, dos=0)

  observeEvent(input$buttonuno,

               {
                 values$uno<-1
                 values$dos<-0
               })

  observeEvent(input$buttondos,

               {
                 values$uno<-0
                 values$dos<-1
               })

  ###################################### View Files Ends Here ########################################





  # Creating pdf and exporting output to pdf

  output$Chemical_Properties <- renderTable({
    if(values$uno || values$dos){

      sink("www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt")
      print(paste("*****************************************************"))
      print(paste("User Information"))
      print(paste("*****************************************************"))
      proj<-input$projID
      name<- input$name
      email<- input$email
      model<- input$model


      print(paste("Project ID:",proj))
      print(paste("Name of investigator:",name))
      print(paste("Name of investigator:",email))
      print(paste("Model Used:",model))
      print(paste("################# End of User Information ##############"))
      sink()
    }

    else
      return()
  })

  # Header for Molecular Formula and Molecular Weight
  output$header_mf_mw<- renderTable({
    if(values$uno || values$dos){
      # Writing Chemical Property
      if(input$Drug_Chemical_Properties)
      {
        header_mf_mw<- "#### Molecular Formula and Molecular Weight ####"
        write.table(header_mf_mw, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })

  ##################################
  # Drug Chemical Properties Tab
  ##################################

  output$chemical_attributes<- renderTable({
    if(values$uno || values$dos){
      # Writing Chemical Property
      if(input$Drug_Chemical_Properties)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        data_chemical_property<-data.frame(MF=MF(read.SDFset(inFile$datapath)[[1]]), MW=MW(read.SDFset(inFile$datapath)[[1]]))

        write.table(data_chemical_property, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })


  ##################################
  # Drug Atomic Properties Tab
  ##################################

  #Header for Atomic attribute
  output$header_atomic_attributes<- renderTable({
    if(values$uno || values$dos){
      if(input$Drug_Atomic_Properties)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        header_atomic_attributes<- "#### Atomic Block ####"
        write.table(header_atomic_attributes, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })

  output$atomic_attributes<- renderTable({
    if(values$uno || values$dos){
      if(input$Drug_Atomic_Properties)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        data_atomic_property<- atomblock(read.SDFset(inFile$datapath)[[1]])
        write.table(data_atomic_property, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })

  #Header for Bond attribute
  output$header_bond_attributes<- renderTable({
    if(values$uno || values$dos){
      if(input$Drug_Atomic_Properties)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        header_bond_attributes<- "#### Bond Block ####"
        write.table(header_bond_attributes, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })

  output$bond_attributes<- renderTable({
    if(values$uno || values$dos){
      if(input$Drug_Atomic_Properties)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        bond_property<- bondblock(read.SDFset(inFile$datapath)[[1]])
        write.table(bond_property, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })




  ##################################
  # Bond Attributes Tab
  ##################################

  #Header for Ring Number
  output$header_ring_attributes<- renderTable({
    if(values$uno || values$dos){
      if(input$Bond_Attributes)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        header_ring_attributes<- "#### Total number of rings  ####"
        write.table(header_ring_attributes, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })


  #Total Ring Attribute
  output$total_ring_number<- renderTable({
    if(values$uno || values$dos){
      if(input$Bond_Attributes)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        total_ring_number<-   rings(read.SDFset(inFile$datapath)[[1]],type="count",upper=1000, arom=TRUE)
        write.table(total_ring_number, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = TRUE, append = T)
      }
      else{
        return()
      }
    }
  })

  #Header for Ring Type
  output$header_ring_type<- renderTable({
    if(values$uno || values$dos){
      if(input$Bond_Attributes)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        header_ring_type<- "#### Ring Aromaticity  ####"
        write.table(header_ring_type, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })


  #Ring Aromaticity
  output$ring_aromaticity<- renderTable({
    if(values$uno || values$dos){
      if(input$Bond_Attributes)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        ring_type<-    rings(read.SDFset(inFile$datapath)[[1]], upper=Inf, type="all", arom=TRUE, inner=FALSE)
        write.table(ring_type[2], file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = TRUE, append = T)
      }
      else{
        return()
      }
    }
  })


  #Header for Ring Structure
  output$header_ring_structure<- renderTable({
    if(values$uno || values$dos){
      if(input$Bond_Attributes)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        header_ring_structure<- "#### Ring Structure  ####"
        write.table(header_ring_structure, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })


  #Ring Structure
  output$ring_structure<- renderTable({
    if(values$uno || values$dos){
      if(input$Bond_Attributes)
      {
        inFile <- input$file
        if (is.null(inFile))
          return(NULL)

        rc <- rings(read.SDFset(inFile$datapath)[[1]], upper=Inf, type="all", arom=TRUE, inner=FALSE)
        rc_df<- c(rc$RINGS[seq(from= 1, to = length(rc$RINGS))])
        rc_df <- as.character(rc_df)
        rc_df[rc_df==""] <- "NA"
        rc_df <- as.factor(rc_df)
        write.table(rc_df, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = TRUE, append = T)
      }
      else{
        return()
      }
    }
  })

  #Header for Bond Matrix
  output$header_bond_matrix<- renderTable({
    if(values$uno || values$dos){
      if(input$Bond_Attributes)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        header_bond_matrix<- "#### Bond Matrix  ####"
        write.table(header_bond_matrix, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })


  #Bond Matrix
  output$bond_matrix<- renderTable({
    if(values$uno || values$dos){
      if(input$Bond_Attributes)
      {
        inFile <- input$file
        if (is.null(inFile))
          return(NULL)

        bond_matrix<- conMA(read.SDFset(inFile$datapath)[[1]], exclude=c("H"))
        write.table(bond_matrix, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = TRUE, append = T)
      }
      else{
        return()
      }
    }
  })

  ##################################
  # Drug Classification Tab
  ##################################
  #Header for Drug Classification
  output$header_drug_classification<- renderTable({
    if(values$uno || values$dos){
      if(input$Drug_Classification)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        header_drug_classification<- "#### Drug Classification  ####"
        write.table(header_drug_classification, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })


  # Drug Classification for output
  output$file_output_drug_classification = renderTable({
    if(values$uno || values$dos){
      if(input$Drug_Classification)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        sdf_input<-read.SDFset(inFile$datapath)[[1]]
        write.SDF(sdf_input,"input.sdf")
        dc<-drug.class.deterministic("input.sdf")

        drug.class.output<-print(paste("The drug is likely to belong", dc, "class."))
          write.table(drug.class.output, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)

      }}
  }, rownames = FALSE, colnames = FALSE,
  include.rownames=FALSE)


  ##################################
  # Drug Neighbors Tab
  ##################################
  #Header for Drug Classification
  output$header_drug_neighbors<- renderTable({
    if(values$uno || values$dos){
      if(input$Nearest_drugs)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        header_drug_neighbors<- "#### Drug Nearest Neighbors  ####"
        write.table(header_drug_neighbors, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })


  #Drug Nearest Neighbor output to the file

  output$file_output_nearest_neighbors = renderTable({
    if(values$uno || values$dos){
      if(input$Nearest_drugs)
      {
        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        sdf_input<-read.SDFset(inFile$datapath)[[1]]
        write.SDF(sdf_input,"input.sdf")
        dc<-drug.class.deterministic("input.sdf")

        organism <- input$organism
        if(organism == "Escherichia coli")
        {
        db<- read.csv2("db/Antibiotic_data_set_Ecoli.csv", header = TRUE, sep = ",")
        resistome_deterministic = data.frame(subset(db, Drug_class == dc))
        drug.neighbor<-data.frame(resistome_deterministic[,c(1,6)])
        nearest.neighbor<-subset(drug.neighbor, !duplicated(Drug_Name))
        write.table(nearest.neighbor, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
        }
        else
        {
          db<- read.csv2("db/Antibiotic_data_set_Paeruginosa.csv", header = TRUE, sep = ",")
          resistome_deterministic = data.frame(subset(db, Drug_class == dc))
          drug.neighbor<-data.frame(resistome_deterministic[,c(1,6)])
          nearest.neighbor<-subset(drug.neighbor, !duplicated(Drug_Name))
          write.table(nearest.neighbor, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
        }
        }}
  },
  include.rownames=FALSE)


  #####################################
  # Resistance Genes Tab
  #####################################
  #Header for Resistance Gene List
  output$header_resistance_genes<- renderTable({
    if(values$uno || values$dos){
      if(input$Resistance_Gene_List)
      {

        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        header_resistance_genes<- "#### Resistance Causing Genes  ####"
        write.table(header_resistance_genes, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
      }
      else{
        return()
      }
    }
  })

  #Resistance_genes

  output$file_output_resistance_genes = renderTable({
    if(values$uno || values$dos){
      if(input$Resistance_Gene_List)
      {
        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        sdf_input<-read.SDFset(inFile$datapath)[[1]]
        write.SDF(sdf_input,"input.sdf")
        dc<-drug.class.deterministic("input.sdf")

        organism <- input$organism
        if(organism == "Escherichia coli")
        {
        db<- read.csv2("db/Antibiotic_data_set_Ecoli.csv", header = TRUE, sep = ",")
        resistome_deterministic = data.frame(subset(db, Drug_class == dc))
        nearest.neighbor<-subset(resistome_deterministic, !duplicated(Drug_Name))
        write.table(nearest.neighbor, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
        }else
        {
          db<- read.csv2("db/Antibiotic_data_set_Paeruginosa.csv", header = TRUE, sep = ",")
          resistome_deterministic = data.frame(subset(db, Drug_class == dc))
          nearest.neighbor<-subset(resistome_deterministic, !duplicated(Drug_Name))
          write.table(nearest.neighbor, file = 'www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt', row.names = FALSE, append = T)
        }

      }}
  },
  include.rownames=FALSE)



  # Printing output file
  output$Properties_output <- renderUI({
    if(values$uno){
      PDFfile="uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt"
      tags$iframe(title = "uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output",
                  src=PDFfile,
                  width="80%",
                  height="400px")
    }
    else
      return()
  })

  # Download printing

  output$downloadData <- downloadHandler(
    filename <- function() {
      paste("uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output", ".txt", sep=".")
    },

    content <- function(file) {
      file.copy("www/uCARE.Chem.Suit.Version.0.2.Chemical.Properties.Output.txt", file)
    },
    contentType = "application/text"
  )



  # Printing Chemical structure plot
  output$Chemical_structure_plot <- renderUI({
    if(values$uno || values$dos){

      if(input$Drug_Chemical_Properties)
      {
        pdf("www/uCARE.Chem.Suit.Version.0.2.Chemical.Structure.Output.pdf")
        inFile <- input$file
        plot(read.SDFset(inFile$datapath)[[1]], print=FALSE, atomnum=TRUE, no_print_atoms="H", main="Chemical Structure")

        graphics.off()


        PDFfile="uCARE.Chem.Suit.Version.0.2.Chemical.Structure.Output.pdf"
        tags$iframe(title = "uCARE.Chem.Chemical.Structure.Output",
                    src=PDFfile,
                    width="80%",
                    height="400px")

      }}

    else
      return()
  })

  ###############################
  # Complete Database

  output$db<- renderDataTable({
    organism <- input$organism
    if(organism == "Escherichia coli")
    {
      db<- read.csv2("db/Antibiotic_data_set_Ecoli.csv", header = TRUE, sep = ",")
    }else
    {
      db<- read.csv2("db/Antibiotic_data_set_Paeruginosa.csv", header = TRUE, sep = ",")
    }


    pubchem <- paste0('<a href="https://pubchem.ncbi.nlm.nih.gov/compound/',db$Drug_Pubchem_ID,'" target="_blank">',db$Drug_Pubchem_ID,'</a>')

    x <- unlist(strsplit(as.character(db$Bibliography),":"))
    x <- trimws(x[seq(2,length(x),2)])

    db_main <- cbind(db[,-c(3,6)], "Drug_Pubchem_ID" = pubchem,"Bibliography"=paste0('<a href="https://www.ncbi.nlm.nih.gov/pubmed/',x,'" target="_blank">',x,'</a>') )


  },escape=F)

  # Complete database finishes here

  ##############################################################
  ####################### Stochastic Tabset Panel Starts here
  ##############################################################
  #Stochastic drug classification

  output$scclass = renderTable({

    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    sdf_input<-read.SDFset(inFile$datapath)[[1]]
    write.SDF(sdf_input,"input.sdf")

      drug.class.sto<- drug.class.stochastic("input.sdf", input$neighbor, input$threshold)
      print(paste("The drug is likely to belong", drug.class.sto[2],"class"))
  },escape = FALSE, include.rownames=FALSE, include.colnames=FALSE)

  # Extracting Drug class
  output$ap_search <- renderDataTable({

    # Reading the Database
    antibiotics <- read.SDFset("db/all_sdf_names.sdf")
    cid(antibiotics) <-  makeUnique(sdfid(antibiotics))
    apset <- sdf2ap(antibiotics)

    # Reading Query
    inFile <- input$file
    sdf_input <- read.SDFset(inFile$datapath)[[1]]
    sdf_input_ap<- sdf2ap(sdf_input)



    df<-cmp.search(apset, sdf_input_ap, type=3, cutoff = 76, quiet=TRUE)

    colnames(df)<- c("Index","Drug_name","Score")
    df

  })


  ##############################################################
  ####################### Stochastic Tabset Panel Starts here
  ##############################################################

  model<- reactive({
    radioButtons("feature","Feature",list("Atom Pair", "Fingerprint"), inline = TRUE)
  })

  clus_type<- reactive({
    radioButtons("clus_type","Plot Type",list("Tree", "Circle"), inline = TRUE)
  })



  output$clus <- renderUI({

    # Reading the database
    antibiotics <- read.SDFset("db/all_sdf_names.sdf")

    # Reading the query
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    sdf_input<-read.SDFset(inFile$datapath)[[1]]


    #Concatanate query and db objects
    antibiotics[[78]]<- sdf_input

    unique_ids<-capture.output(cid(antibiotics) <-  makeUnique(sdfid(antibiotics)))
    no_duplicate_compound<-length(grep("appended", unique_ids, ignore.case = TRUE))
    if(no_duplicate_compound==0)
    {
      apset <- sdf2ap(antibiotics)

      apdups <- cmp.duplicated(apset, type=1)
      antibiotics[which(!apdups)]; apset[which(!apdups)]

      if(input$feature == "Atom Pair" && input$clus_type== "Tree")
      {
        dummy <- cmp.cluster(db=apset, cutoff=0, save.distances="distmat.rda", quiet=TRUE)
        load("distmat.rda")
        hc <- hclust(as.dist(distmat), method="single")
        hc[["labels"]] <- cid(apset) # Assign correct item labels

        pdf(file = "www/query_db_cluster.pdf")
        par(mar=c(7,5,1,5),cex=0.5, cex.main=1.5)
        plot(as.dendrogram(hc), edgePar=list(col=1, lwd=2), horiz=F)
        title("Binning clustering with atom pairs as features (Tree)", line = -3)
        rect.hclust(hc, h=mean(hc[[2]]), border=c(2:8))
        dev.off()
      }

      else if(input$feature == "Atom Pair" && input$clus_type== "Circle")
      {
        dummy <- cmp.cluster(db=apset, cutoff=0, save.distances="distmat.rda", quiet=TRUE)
        load("distmat.rda")
        hc <- hclust(as.dist(distmat), method="single")
        hc[["labels"]] <- cid(apset) # Assign correct item labels
        pdf(file = "www/query_db_cluster.pdf")

        colors = c(2:30)
        clus4 = cutree(hc, h=mean(hc[[2]]))
        plot(as.phylo(hc), type = "fan", tip.color = colors[clus4],
             #label.offset = 1,
             cex = 0.7)
        title("Binning clustering with atom pairs as features (Circle)", line = 1)
        dev.off()
      }


      else if(input$feature == "Fingerprint" && input$clus_type== "Tree")
      {
        # Reading the database
        antibiotics_db <- read.SDFset("db/all_sdf_names.sdf")

        # Reading the query
        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        sdf_input<-read.SDFset(inFile$datapath)[[1]]

        write.SDF(sdf_input,"input.sdf")

        # Test run
        sdf_input_fp<-read.SDFset(sdfstr = "input.sdf")

        unique_ids<-capture.output(cid(antibiotics_db) <-  makeUnique(sdfid(antibiotics_db)))
        cid(sdf_input_fp)<- sdfid(sdf_input_fp)

        apset_db <- sdf2ap(antibiotics_db)
        fpset_db <- desc2fp(apset_db)

        apset_in <- sdf2ap(sdf_input_fp)
        fpset_in <- desc2fp(apset_in)

        concatenated_fpset<-c(fpset_db[1:77], fpset_in[1])

        simMA <- sapply(cid(concatenated_fpset), function(x) fpSim(concatenated_fpset[x], concatenated_fpset, sorted=FALSE))
        hc <- hclust(as.dist(1-simMA), method="single")

        pdf(file = "www/query_db_cluster.pdf")
        par(mar=c(7,5,1,5), cex.main=1, cex=0.7)
        plot(as.dendrogram(hc), edgePar=list(col=1, lwd=2), horiz=F)
        title("Binning clustering with atom fingerprints as features (Tree)", line = -3)
        rect.hclust(hc, h=mean(hc[[2]]), border=c(2:8))
        dev.off()

      }

      else if(input$feature == "Fingerprint" && input$clus_type== "Circle")
      {
        # Reading the database
        antibiotics_db <- read.SDFset("db/all_sdf_names.sdf")

        # Reading the query
        inFile <- input$file
        if (is.null(inFile))
          return(NULL)
        sdf_input<-read.SDFset(inFile$datapath)[[1]]

        write.SDF(sdf_input,"input.sdf")

        # Test run
        sdf_input_fp<-read.SDFset(sdfstr = "input.sdf")


        unique_ids<-capture.output(cid(antibiotics_db) <-  makeUnique(sdfid(antibiotics_db)))
        cid(sdf_input_fp)<- sdfid(sdf_input_fp)

        apset_db <- sdf2ap(antibiotics_db)
        fpset_db <- desc2fp(apset_db)

        apset_in <- sdf2ap(sdf_input_fp)
        fpset_in <- desc2fp(apset_in)

        concatenated_fpset<-c(fpset_db[1:56], fpset_in[1])

        simMA <- sapply(cid(concatenated_fpset), function(x) fpSim(concatenated_fpset[x], concatenated_fpset, sorted=FALSE))
        hc <- hclust(as.dist(1-simMA), method="single")

        pdf(file = "www/query_db_cluster.pdf")

        colors = c(2:30)
        clus4 = cutree(hc, h=mean(hc[[2]]))
        plot(as.phylo(hc), type = "fan", tip.color = colors[clus4],
             #label.offset = 1,
             cex = 0.7)
        title("Binning clustering with atom fingerprints as features (Circle)", line = 1)


        dev.off()

      }

    }

    PDFfile="query_db_cluster.pdf"
    tags$iframe(title = "query_db_cluster",
                src=PDFfile,
                width="90%",
                height="500px")

  })

  #### Counter #####
  output$counter<- renderTable({
    counter<- read.table("www/counter.csv", header = FALSE)
    present_count<- counter$V1
    added_count<- present_count + 1
    final_count<-added_count
    write.csv(added_count,file='www/counter.csv')
    tail(final_count,n=1)

  },include.colnames=FALSE)




  # FAQ
  output$faq<- renderUI({

    tags$ol(br(),
            h2("Frequently Asked Question (FAQ)"),br(),

            # FAQ header
            h4("#Introduction to uCARE Chem Suite Version 0.2"),

            # Question 1
            tags$a(href = "#A", "Q 1. What is uCARE Chem Suite Version 0.2?"),
            tags$a(name = "#A"), h5("A 1. uCARE Chem Suite Version 0.2 is an online chemoinformatics suite written in R.
                                    It is based on homology based approach to predict bacterial resistome profile for input
                                    chemical compound structures"),
            # Question 2
            tags$a(href = "#B", "Q 2. What is the rationale behind developing uCARE Chem Suite Version 0.2?"),
            tags$a(name = "#B"), h5("Functional annotation of proteins on the basis of homology is a common practice
                                (Loewenstein", tags$i("et al.,")," 2009; Mahlich", tags$i("et al.,"),"  2018).
                                Moreover, Soares and Carneigo, (2002) suggested that due to the similarity of
                                drugs", tags$i("i.e.,"), " the similarity in terms of similarity of chemical structure,
                                the similarity in the mechanism of action, and similarity of pharmaceutical effect, similar
                                drugs are interchangeably used (The drug class effect).",br(),

                                "Therefore, uCARE Chem Suite predicts the bacterial resistome on the basis of
                                the hypothesis that, the drug class effect is not limited to pharmaceutical properties but
                                will also encompass the resistome profile", tags$i("i.e.,"),  "resistome profile prediction of unknown
                                small molecules on the basis of its structural homology with antibiotics of known resistome."),
            # Question 3
            tags$a(href = "#C", "Q 3. What are different models offered by uCARE Chem Suite Version 0.2 for similarity search?"),br(),
            tags$a(name = "#C"), h5("A 3. Presently, the suite utilizes two models", tags$i("viz."), "deterministic model based upon literature driven rules
                                and k- Nearest Neighbor based stochastic model."),
            #Question 4
            tags$a(href = "#D", "Q 4. What is uCARE Chem Suite Version 0.2's deterministic model?"),br(),
            tags$a(name = "#D"), h5("A 4. uCARE Chem Suite Version 0.2's deterministic model determines the class of compound
                                    on the basis of algorithm formulated by exhaustive manual annotation of resistant drug classes."),

                                tags$i(p("Example. 1. Any compound is likely to be a beta lactam drug if it has a beta lactam ring i.e.
                                if the drug compound contains an alipahtic ring consiting of 3 carbon atoms and one nitrogen atom, it is
                                likely to be a Beta Lactam.")),
            #Question 5
            tags$a(href = "#E", "Q 5. What is uCARE Chem Suite Version 0.2's Stochastic model?"),br(),
            tags$a(name = "#E"), h5("A 5. uCARE Chem Suite Version 0.2's Stochastic model provides information of similar compounds by aligning
                                    query compound to the database compounds using tanimoto coeffiecient of atom pairs as a similarity measure."),

            br(),h4("#Usage prerequisite"),
            #Question 6
            tags$a(href = "#F", "Q 6. What are different structural formats supported by uCARE Chem Suit Version 0.2? "),br(),
            tags$a(name = "#F"), h5("A 6. uCARE Chem Suite Version 0.2 presently works on structure data files (SDF) only."),

            br(),h4("#Usage"),
            #Question 7
            tags$a(href = "#G", "Q 7. What are different outputs that one can get using uCARE Chem Suite Version 0.2's deterministic model?"),br(),
            tags$a(name = "#G"), h5("A 7. uCARE Chem Suite Version 0.2's deterministic model provide utilities to visualize molecular
                                    structure and formula of the compound (Chemical Property tab); atom blocks and compound blocks of the compound
                                    (Atomic Property tab); ring number, ring type, ring structure and ring coordinates of the compound (Bond Attributes tab);
                                    classification of compound among predefined classes on the basis of its physiochemical property (Drug Classification tab);
                                    other compounds of the same class from the database (Nearest Drug/s tab); list of genes reported to confer resistance to similar
                                    compounds (Resistance Gene list tab) and Export/Download tab to export the output in pdf and text format."),
            #Question 8
            tags$a(href = "#H", "Q 8. What are different outputs that one can get using uCARE Chem Suite Version 0.2's Stochastic model?"),br(),
            tags$a(name = "#H"), h5("A 8. uCARE Chem Suit Version 0.2's Stochastic model predicts the class of input compound (Drug classification tab)
                                    using k-nearest neighbor algorithm, and provides utilities to align query
                                    structure against the database compounds (Database Query Search tab); perform clustering of query chemical structure
                                    against the database compounds by binning algorithm using atom pair and fingerprints as feature (Query-DB clustering tab)
                                    and access the entire database through ease to use graphical user interface (Database tab).
                                    ", br(),br(),
                                    "Note: The query compound is labelled as CID number (Number) in the clustering plots."),
            #Question 9
            tags$a(href = "#I", "Q 9. I am in Output console; how can I go to the HOME console?"),br(),
            tags$a(name = "#I"), h5("A 9. Just referesh the tab!"),

            br(),h4("#Errors, bug fix and feedback"),

            #Question 10
            tags$a(href = "#J", "Q 10. I am submitting second molecule in same session and getting error in Export tab's Resistance Gene list viz. 'Drug_Name' object not found. What should I do?"),br(),
            tags$a(name = "#J"), h5("A 10. Presently, our server works on single instance. Just referesh the tab!"),

            #Question 11
            tags$a(href = "#K", "Q 11. Query-DB clustering tab showing Error 404. What should I do? "),br(),
            tags$a(name = "#K"), h5("A 11. The server might be busy, so please try later after some time."),

            #Question 12
            tags$a(href = "#L", "Q 12. Is it possible to see the main algorithms running behind this server?"),br(),
            tags$a(name = "#L"), h5("A 12. You are most welcome to download the entire code! The scripts of prediction algorithm
                                    can be browsed from:", a("GitHub Repository", href="https://github.com/sauravbsaha/uCAREChemSuiteCLI/")),

            #Question 13
            tags$a(href = "#M", "Q 13. I have found a bug or/and want to provide feedback. Where should I report it?"),br(),
            tags$a(name = "#M"), h5("A 13. Please shoot a mail to saurav.saha@shiats.edu.in. Your feedbacks are most welcome.")
        )
  })





  #########################################################################
  ################### Rule based analysis Ends Here #################################
  #########################################################################
  about<- helpText(br(),
                   p("uCARE Chem Suite is a R package to predict the
                     resistome profile of ", tags$i("Escherichia coli (E. coli)"),"and",
                      tags$i("Pseudomonas aeruginosa (P. aeruginosa)"),  "towards
                     the novel candidate drug. In addition to resistome prediction,
                     the suite enables the user to visualize the chemical structure,
                     classify compound in predefined 19 drug classes, perform pairwise
                     alignment and cluster with database compounds using graphical user
                     interface."),

                   p("The usage can be browsed through FAQ page."),

                   HTML('<center><img src="drug_wordcloud4.png" width="750"></center>')



                   )
  db_query_Search<- helpText(
    h3("Similarity search result"),
    p("Similarity scores are aligned in descending order. Higher the score, higher is the query subject similarity.",br(),
      "Note: Class annotations of the nearest neighbor can be found in database tab")


  )

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Main Panel Starts Here $$$$$$$$$$$$$$$$$$$$$$$
  #If-else statement to print on main panel

  output$home <- renderUI({
    if(is.null(input$file))
      tabsetPanel(tabPanel("Home",about,br(),h4("Usage Counter:")
                           ,tableOutput("counter")
      ),

      tabPanel("FAQ", uiOutput("faq"),br()),
      tabPanel("Contact Us", h2("Team"),tags$img(src='Contact.jpg', height=300, width=750) )
      )



    else
    {
      if(input$application =="Visualization")
      {
        tabsetPanel(tabPanel("Chemical Properties", h3("Chemical Properties"), tableOutput("mf"), h3("Chemical structure"), plotOutput("img", width = "550px", height = "500px")),

                    tabPanel("Atomic Properties", h3("Atom block"), tableOutput("mytable"),
                             h3("Bond Block"), tableOutput("bondblock")),

                    tabPanel("Bond Attributes", h3("Total number of Rings"), tableOutput("ring"),h3("Ring Types"), tableOutput("rnum"), h3("Ring structures"), tableOutput("rstructure"), h3("Bond matrix"), tableOutput("bmatrix"))
        )

      }
      else if(input$model == "Deterministic Model" && input$application == "Resistome Prediction")
      {
        tabsetPanel(tabPanel("Drug Classification", h3("Classification"),h4(tableOutput("cclass")) ),

                    tabPanel("Nearest drug/s", h3( textOutput("mymodel")),h4("Most Nearest drug accroding to Drug structure"), tableOutput("dclass") ),

                    tabPanel("Resistance Gene List",h3("Predicted resistant gene/s"),  dataTableOutput("gclass")),

                    tabPanel("Export/Downloads", br(),h3("Select output from the following"),
                             checkbox_Drug_Chemical(),checkbox_Drug_Atomic(),checkbox_Bond_Attributes(),
                             checkbox_Drug_Classification(),checkbox_Nearest_drug(), checkbox_Resistance_Gene_List(),
                             export_submit(),
                             h3("Generating output file"),"First click Create followed by Download button" ,
                             br(),br(),export_download() , downloadButton("downloadData", label = "Download"),br(),br(),

                             # h3("Output File"),


                             # Creating file
                             tableOutput("Chemical_Properties"),

                             ##################################
                             # Drug Chemical Properties Tab
                             ##################################

                             #Header for Chemical Attributes
                             tableOutput("header_mf_mw"),

                             #Printing Chemical Attributes
                             tableOutput("chemical_attributes"),

                             ##################################
                             # Drug Atomic Properties Tab
                             ##################################

                             #Header for Atomic Attributes
                             tableOutput("header_atomic_attributes"),

                             #Printing Atomic Attributes
                             tableOutput("atomic_attributes"),

                             #Header for Bond Attributes
                             tableOutput("header_bond_attributes"),

                             #Printing Bond Attributes
                             tableOutput("bond_attributes"),

                             ##################################
                             # Bond Attributes Tab
                             ##################################
                             #Header Bond Attributes
                             tableOutput("header_ring_attributes"),

                             # Printing Total ring number
                             tableOutput("total_ring_number"),

                             #Header Ring Type
                             tableOutput("header_ring_type"),

                             #Ring Type
                             tableOutput("ring_aromaticity"),

                             #Header Ring Structure
                             tableOutput("header_ring_structure"),

                             #Ring Structure
                             tableOutput("ring_structure"),

                             #Header Bond Matrix
                             tableOutput("header_bond_matrix"),

                             #Ring Matrix
                             tableOutput("bond_matrix"),

                             ##################################
                             # Drug Classification Tab
                             ##################################
                             #Header Drug Calssification Tab
                             tableOutput("header_drug_classification"),

                             #Drug Calssification Tab
                             tableOutput("file_output_drug_classification"),

                             ##################################
                             # Neighboring Drugs Tab
                             ##################################
                             #Header Drug Neighbors Tab
                             tableOutput("header_drug_neighbors"),

                             #Drug Nearest Neighbors
                             tableOutput("file_output_nearest_neighbors"),

                             ##################################
                             # Resistance Genes Tab
                             ##################################

                             #Header Drug Neighbors Tab
                             tableOutput("header_resistance_genes"),

                             #Resistance Gene List
                             tableOutput("file_output_resistance_genes"),

                             ##################################
                             # Chemical Structure
                             ##################################
                             #Printing Text output
                             uiOutput("Properties_output"),


                             #Printing Plot output
                             #h3("Chemical Structure"),
                             uiOutput("Chemical_structure_plot"),
                             br() )

        )}
      else if (input$model == "Stochastic Model" && input$application == "Resistome Prediction"){
        tabsetPanel(
          tabPanel("Drug Classification", h3("Classification"),h4(tableOutput("scclass"))),
          tabPanel("Database Query Search",br(), db_query_Search,br(), dataTableOutput("ap_search")),

          tabPanel("Query-DB Clustering",
                   br(),model(),
                   clus_type(),
                   uiOutput("clus"),br(),br(),br()),


          tabPanel("Database", dataTableOutput("db"))
        )

      }

    }
  })


}

#$$$$$$$$$$$$$$$$$$$ Program ends here $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
shinyApp(ui = ui, server = server)
