nextflow.enable.dsl=2

if (!params.obj_key_samplesheet) {
  exit 1, '--obj_key_samplesheet flag empty! Samplesheet containing S3 objects keys not found. Key names should consist of an objects path after the S3 bucket name'
}

if (!params.s3_bucket) {
  exit 1, '--s3_bucket flag empty! Please provide the S3 bucket containing the objects to be tagged'
}

if (!params.clipath) {
  params.clipath = "aws"
}

if (!params.aws_region) {
  params.aws_region = "us-east-1"
}

params.tag_keys = ""
params.tag_values = ""

def make_tagJson(String keys = "", String values = "") {
  def keyList   = keys.split(',').toList()
  def valueList = values.split(',').toList()
    
  if( keyList.size() < valueList.size() ) {
    throw new IllegalArgumentException("Number of key elements must be the same or larger than the number of value elements. " + 
                                        "Key #: ${keyList.size()} vs Value #: ${valueList.size()}")
  }
  
  def tagList = keyList.withIndex().collect { k, i -> [Key: k, Value: ( valueList[i] ?: "")] }
  def tagSet  = [ TagSet: tagList ]
  def tagJson = groovy.json.JsonOutput.toJson(tagSet)
  
  return tagJson
}

// Tag S3 objects using AWS CLI
process put_object_tag {
  container "${params.use_container}"

  input:
  tuple val(obj_key), val(tagJson)

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def obj_key_string = obj_key.toString()

  """
  ${params.clipath} s3api put-object-tagging \\
    --region '${params.aws_region}' \\
    --bucket '${params.s3_bucket}' \\
    --key '${obj_key_string}' \\
    --tagging '${tagJson}'
  """
}

// Get tags of S3 objects using AWS CLI
process get_object_tag {
  input:
  val(obj_key)

  output:
  path("individual_obj_tag_${task.index}.txt")

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def obj_key_string = obj_key.toString()

  """
  ${params.clipath} s3api get-object-tagging \\
    --region '${params.aws_region}' \\
    --bucket '${params.s3_bucket}' \\
    --key '${obj_key_string}' \\
    --output text | awk '{print "'s3://${params.s3_bucket}/${obj_key_string}'\t" \$0}' > "individual_obj_tag_${task.index}.txt"
  """
}

// Concatenate outputs from individual object tags into a single file
process concatenate_outputs {
    input:
    path('individual_obj_tag_*.txt')

    output:
    path('all_obj_tags.txt')

    publishDir '.', mode: 'copy'

    script:
    """
    cat individual_obj_tag_*.txt > all_obj_tags.txt
    """
}

workflow {
  obj_keys_ch = Channel.fromPath(file(params.obj_key_samplesheet))
    .splitCsv()
    .flatten()

  map_obj_with_tag_json_ch = obj_keys_ch
    .map { it -> 
      tuple(it, make_tagJson(params.tag_keys, params.tag_values))
    }

  if (params.tag_keys) {
    map_obj_with_tag_json_ch | put_object_tag
  } else {
    obj_keys_ch | get_object_tag | collect | concatenate_outputs
  }
}
