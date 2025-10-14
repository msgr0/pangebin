include { MAIN               } from './'


def pty = [0, 1]
def tests = ["pan"]



workflow {
    // prendi una cartella "dataset" e itera su tutti i sample in samples.csv, in input

    if (!params.samples) {
        print "Samples csv not provided, trying to load whole dataset"
        samples_ch = Channel.empty()
        dataset_ch = Channel.fromPath(params.dataset)

        samples_ch = dataset_ch.map { d ->
            def dir = file(d)
            if (!dir.isDirectory()) {
                error "Dataset path ${d} is not a directory"
            }
            def samples = dir.list().findAll { it -> file("${dir}/${it}/unicycler.gfa.gz").exists() && file("${dir}/${it}/skesa.gfa.gz").exists() }

            if (samples.size() == 0) {
                error "No samples found in dataset directory ${d}"
            }
            return samples
        }
        samples_ch = samples_ch.flatMap{ it -> it }
        samples_ch = samples_ch.map { sample ->
            path = file("${params.dataset}/${sample}")
            [sample, path]
        }
    }

    else {
        print "Samples csv provided, loading samples from ${params.samples}"
    
        params.samples = file(params.samples)
        samples_ch = Channel.fromPath(params.samples).splitCsv(header:false)
        dataset_ch = Channel.fromPath(params.dataset)
        samples_ch = samples_ch.map { sample ->
                sample = sample[0].toString().trim()
                path = file("${params.dataset}/${sample}")
                [sample, path]
            }
    }


    samples_ch.iew()

    input_ch = samples_ch.map { sample, path ->
        def fmeta = [:];
        fmeta.id = sample;
        def f1 = file("${path}/unicycler.gfa.gz")
        def f2 = file("${path}/skesa.gfa.gz")
        if (f1.exists() && f2.exists()) {
            [[id: sample], [f2, f1]]  // always s first,
        } else {
            error "Need both assemblies for sample ${sample}"
        }
        [fmeta, [f2, f1]]
    }


    MAIN( input_ch)
}
