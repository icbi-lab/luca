// Import generic module functions
include { dump_params_yml; indent_code_block } from "./parametrize"

params.parametrize = true
params.implicit_params = true
params.meta_params = true

process JUPYTERNOTEBOOK {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(notebook)
    val(parameters)
    path(input_files)

    output:
    tuple val(meta), path("*.html"), emit: report
    path("artifacts/*"), emit: artifacts, optional: true
    path "versions.yml"          , emit: version

    script:
    def prefix   = meta.id
    def kernel   = task.ext.kernel ?: '-'

    // Dump parameters to yaml file.
    // Using a yaml file over using the CLI params because
    //  * no issue with escaping
    //  * allows to pass nested maps instead of just single values
    def params_cmd = ""
    def render_cmd = ""
    if (params.parametrize) {
        nb_params = [:]
        if (params.implicit_params) {
            nb_params["cpus"] = task.cpus
            nb_params["artifact_dir"] = "artifacts"
            nb_params["input_dir"] = "./"
        }
        if (params.meta_params) {
            nb_params["meta"] = meta
        }
        nb_params += parameters
        params_cmd = dump_params_yml(nb_params)
        render_cmd = "papermill -f .params.yml"
    } else {
        render_cmd = "papermill"
    }

    """
    set -o pipefail

    # Dump .params.yml heredoc (section will be empty if parametrization is disabled)
    ${indent_code_block(params_cmd, 4)}

    # Create output directory
    mkdir artifacts

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"
    export NUMBA_NUM_THREADS="${task.cpus}"

    # Convert notebook to ipynb using jupytext, execute using papermill, convert using nbconvert
    jupytext --to notebook --output - --set-kernel ${kernel} ${notebook} > ${notebook}.ipynb
    ${render_cmd} ${notebook}.ipynb ${notebook}.executed.ipynb
    jupyter nbconvert --stdin --to html --output ${prefix}.html < ${notebook}.executed.ipynb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        - jupytext: \$( jupytext --version )
        - samtools: \$( python -c "import ipykernel; print(ipykernel.__version__)" )
        - nbconvert: \$(jupyter nbconvert --version)
        - papermill: \$(papermill --version | cut -f1 -d' ')
    END_VERSIONS
    """
}
