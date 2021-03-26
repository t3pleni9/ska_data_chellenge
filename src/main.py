from skavengers import Pipeline, Executor, PipelineSetup

if __name__ == '__main__':
    process_executor = Executor(3)
    pipeline_setupe = PipelineSetup()
    pipeline_config = './skavengers/pipeline_setup.ini'
    pipeline = Pipeline(pipeline_setupe, pipeline_config)

    pipeline.init()

    pipeline.execute(process_executor)
