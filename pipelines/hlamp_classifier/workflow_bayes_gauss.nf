// runs nextflow pipline to find copy number classification from scRNA-seq data
// Reads a yaml files with ALL parameters, runs one config at a time

nextflow.enable.dsl=1

// ------------ Functions ------------

import java.time.LocalDateTime
import java.text.SimpleDateFormat

def randomFileName(rndLen = 10) {
  // Generate a timestamped random file name
  def rndStr = org.apache.commons.lang.RandomStringUtils.random(rndLen, true, true).toString()
  def formatter = new SimpleDateFormat("yyyy-MMM-dd-HH-mm-ss")
  def strDate = formatter.format(new Date())
  return strDate + '-' + rndStr
}

def parse_csv_to_dict_list(infile) {
    // define a function that reads a csv file and converts it to a list of dictionaries
    def header = infile.readLines()[0].split(',')
    content = infile.readLines()[1..-1]
    dictList = content.collect { line ->
        def fields = line.split(',')
        def dict = [:]
        for (int i = 0; i < header.size(); i++) {
            dict[header[i]] = fields[i]
        }
        return dict
    }
    return dictList
}


@Grab('org.yaml:snakeyaml:1.29')
import org.yaml.snakeyaml.Yaml

def parse_yaml_to_dict_list(infile) {
    // Initialize YAML parser
    Yaml yaml = new Yaml()
    // Read the YAML file content
    String content = infile.text
    // Parse the YAML content
    Map parsedYaml = yaml.load(content)
    // Convert the parsed YAML content to a list of dictionaries
    def dictList = []
    parsedYaml.each { key, value ->
        dictList << value
    }
    return dictList
}


def dictToString(dict) {
    // Convert dictionary to a string "--key value" pairs
    // if testing is true, first run dict through set_default_dict
    // Make a copy of the dict
    dict = dict.clone()    
    // remove the id key
    dict.remove('id')
    def str = ""
    for (key in dict.keySet()) {
        str += "--${key} ${dict[key]} "
    }
    return str
}



// ------------ Parameters ------------

// for createPermutedYs
params.configs_path = './batches/local_test/run_configs.yaml'
params.out_dir = '../results/deliverables/'
params.source_dir = '../src/deconfounder'
// Add a test boolean param that will set scresister_max_iters to 100
params.test = false
params.K = 5 // the number of configs to run sequentially on each node

// ------------ Variables ------------
deliverableDir = params.out_dir + workflow.scriptName.replace('.nf','')
// Convert to absolute path 
deliverableDir = (new File(deliverableDir)).getCanonicalPath() + '/' + randomFileName()

// dictList = parse_csv_to_dict_list(file(params.configs_path))
dictList = parse_yaml_to_dict_list(file(params.configs_path))

// if params.test is true, set max_M to 2 subset_counts to 100
if (params.test) {
    // Print a message
    println "Setting max_M to 2 and subset_counts to 100 and n_steps to 100"
    // Change both 
    dictList.each { dict ->
        dict['subset_counts'] = 100
    }
    dictList.each { dict ->
        dict['max_M'] = 2
    }
    dictList.each { dict ->
        dict['n_steps'] = 100
    }
    // Only keep the first 2 configs
    // dictList = dictList[0..1]
    dictList = dictList[0..21]
}


// The output directory, extracted from the dictList
// 
outDirs = dictList.collect { it['outDir'] }
// get the parent direcotry of the parent directory of the first outDir
outDir = new File(outDirs[0]).parentFile.parentFile
summaryOutDir = outDir.getCanonicalPath() + '/' + 'summary'


// ------------ Channels ------------
// Create a channel that goes from 1 to the length of the dictList
sweeps = Channel.from(1..dictList.size())
sweeps.into { sweeps0; sweeps1; sweeps2 }


// ------------ Processes ------------

// Build the code
// Assume being run from with cdm/pipelines
process buildCode {
  cache false
  executor 'local'
  output:
    file 'code' into codeR

  script:
  """
  echo $PWD
  mkdir -p code/
  ln -sf ${params.source_dir}/mixture_model/*.py code
  """
}

codeR.into {
    codeR0
    codeR1
    codeR2
    codeR3
}

// For each config (row of the csv)
// Run the workflow pipeline (call infercnv clones, make some rudimentary plots, and create nextflow script for the next bout)

process ClassifyCNV {
    // errorStrategy { ((task.exitStatus == 1 || task.exitStatus == 140 || task.exitStatus == 137 || task.exitStatus == 130) && task.attempt <= maxRetries) ? 'retry' : 'ignore' }
    memory { 20.GB + 2.GB * task.attempt }
    cpus 1
    maxRetries 0
    time { 0.h + 1.h * task.attempt }
    queue 'componc_cpu'

    input:
        file codeR0
        // Get a chunk (list) of integers from the sweeps channel
        // val chunk from sweeps1.collate(params.K)
        val chunk from sweeps1.collate(params.K, params.K)

    output:
        // Generate a dummy output file for each integer in the chunk
        file('output_*.txt') into outs_Xprime

    // (For simplicity, publish to a common directory)
    publishDir "${deliverableDir}/common", mode: 'copyNoFollow', overwrite: true

    script:
    """
    echo "Processing chunk: ${chunk}"
    # Convert the Groovy list to a Bash array
    CHUNK=(${chunk.join(' ')})
    for id in "\${CHUNK[@]}"; do
       id0=\$(( id - 1 ))
       echo "Processing id=\${id0}"
       if [ "${params.test}" = true ]; then
            python3 code/mixture_models_pyro_gauss_chunk.py ${params.configs_path} \${id0} --test
       else
            python3 code/mixture_models_pyro_gauss_chunk.py ${params.configs_path} \${id0}
       fi
       touch output_\${id0}.txt
    done
    """
}


// Collect all output files into a single list when ClassifyCNV is done
finalOutputs = outs_Xprime.collect()

process FinalProcess {
    errorStrategy { ((task.exitStatus == 1 || task.exitStatus == 140 || task.exitStatus == 137 || task.exitStatus == 130) && task.attempt <= maxRetries) ? 'retry' : 'ignore' }
    memory { 10.GB + 2.GB * task.attempt }
    cpus 1
    maxRetries 0
    time { 0.h + 1.h * task.attempt }
    queue 'cpushort'

    input:
        // This will be a single list containing all output file names
        file codeR1
        val files from finalOutputs
    output:
        file('final.txt') into finalOut

    script:
    """
    # Get the dirname of the dirname of the outDir for the first config file
    # if params.test is true, pass --test to the python script
    python3 code/mixture_utils_files.py ${outDir}  ${summaryOutDir} --model gauss
    touch final.txt
    """
}


// Write pipeline info to a file
process summarizePipeline {
    executor 'local'
    cache false
    output:
        file 'pipeline-info.txt'   
    publishDir deliverableDir, mode: 'move', overwrite: true
    """
    # Get the dirname of the dirname of the outDir for the first config file
    echo 'outDir':  ${outDir} >> pipeline-info.txt
    echo 'SummaryDir' ${summaryOutDir} >> pipeline-info.txt
    echo 'scriptName: $workflow.scriptName' >> pipeline-info.txt
    echo 'start: $workflow.start' >> pipeline-info.txt
    echo 'runName: $workflow.runName' >> pipeline-info.txt
    echo 'nextflow.version: $workflow.nextflow.version' >> pipeline-info.txt
    """
}