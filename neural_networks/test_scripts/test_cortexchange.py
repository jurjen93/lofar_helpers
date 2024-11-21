from cortexchange.wdclient import init_downloader

init_downloader(
    url="https://researchdrive.surfsara.nl/public.php/webdav/",
    login="JeofximLVcr8Ttm",
    password="?CortexAdminTest1?",
    cache="/home/larsve/.cache/cortexchange",
    # cache=".cache/cortexchange",
)

from cortexchange.architecture import get_architecture, Architecture

TransferLearning: type(Architecture) = get_architecture("surf/TransferLearning")
model = TransferLearning(device="cpu", model_name="surf/efficientnet_october_09535")

# torch_tensor = model.prepare_data(
#     "ILTJ160454.72+555949.7_selfcal/selfcal_007-MFS-image.fits"
# )
torch_tensor = model.prepare_data(
    "/scratch-shared/CORTEX/public.spider.surfsara.nl/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/stop/ILTJ142906.77+334820.3_image_009-MFS-image.fits"
)
print(torch_tensor.shape)
result = model.predict(torch_tensor)
print(result)
