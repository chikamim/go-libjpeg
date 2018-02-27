package jpeg

/*
#include <stdio.h>
#include <stdlib.h>
#include "jpeglib.h"
#include "jpegint.h"
#include "jpeg.h"

void error_panic(j_common_ptr cinfo);

static struct jpeg_compress_struct *new_compress(void) {
	struct jpeg_compress_struct *cinfo = (struct jpeg_compress_struct *)calloc(sizeof(struct jpeg_compress_struct), 1);
	struct jpeg_error_mgr *jerr = (struct jpeg_error_mgr *)calloc(sizeof(struct jpeg_error_mgr), 1);

	jpeg_std_error(jerr);
	jerr->error_exit = (void *)error_panic;
	jpeg_create_compress(cinfo);
	cinfo->err = jerr;

printf("HELLO");
fprintf(stderr, "HELLO ERR");

	return cinfo;
}

static void destroy_compress(struct jpeg_compress_struct *cinfo) {
	free(cinfo->err);
	jpeg_destroy_compress(cinfo);
	free(cinfo);
}

static void encode_gray(j_compress_ptr cinfo, JSAMPROW pix, int stride) {
	// Allocate JSAMPIMAGE to hold pointers to one iMCU worth of image data
	// this is a safe overestimate; we use the return value from
	// jpeg_read_raw_data to figure out what is the actual iMCU row count.
	JSAMPROW *rows = alloca(sizeof(JSAMPROW) * ALIGN_SIZE);

	int v = 0;
	for (v = 0; v < cinfo->image_height; ) {
		// First fill in the pointers into the plane data buffers
		int h = 0;
		for (h = 0; h < DCTSIZE * cinfo->comp_info[0].v_samp_factor; h++) {
			rows[h] = &pix[stride * (v + h)];
		}
		// Get the data
		v += jpeg_write_raw_data(cinfo, &rows, DCTSIZE * cinfo->comp_info[0].v_samp_factor);
	}
}

static void encode_rgba(j_compress_ptr cinfo, JSAMPROW pix, int stride) {
	JSAMPROW rows[1];

	int v;
	for (v = 0; v < cinfo->image_height; ) {
		rows[0] = &pix[v * stride];
		v += jpeg_write_scanlines(cinfo, rows, 1);
	}
}


static void encode_ycbcr(j_compress_ptr cinfo, JSAMPROW y_row, JSAMPROW cb_row, JSAMPROW cr_row, int y_stride, int c_stride, int color_v_div) {
	// Allocate JSAMPIMAGE to hold pointers to one iMCU worth of image data
	// this is a safe overestimate; we use the return value from
	// jpeg_read_raw_data to figure out what is the actual iMCU row count.
	JSAMPROW *y_rows = alloca(sizeof(JSAMPROW) * ALIGN_SIZE);
	JSAMPROW *cb_rows = alloca(sizeof(JSAMPROW) * ALIGN_SIZE);
	JSAMPROW *cr_rows = alloca(sizeof(JSAMPROW) * ALIGN_SIZE);
	JSAMPARRAY image[] = { y_rows, cb_rows, cr_rows };

	int v = 0;
	for (v = 0; v < cinfo->image_height; ) {
		int h = 0;
		// First fill in the pointers into the plane data buffers
		for (h = 0; h <  DCTSIZE * cinfo->comp_info[0].v_samp_factor; h++) {
			y_rows[h] = &y_row[y_stride * (v + h)];
		}
		for (h = 0; h <  DCTSIZE * cinfo->comp_info[1].v_samp_factor; h++) {
			cb_rows[h] = &cb_row[c_stride * (v / color_v_div + h)];
			cr_rows[h] = &cr_row[c_stride * (v / color_v_div + h)];
		}
		// Get the data
		v += jpeg_write_raw_data(cinfo, image, DCTSIZE * cinfo->comp_info[0].v_samp_factor);
	}
}

void print_cinfo(j_compress_ptr cinfo) {
  if (!cinfo) {
    fprintf(stderr,"print_cinfo called with NULL cinfo\n");
  }

  fprintf(stderr, "printing cinfo...\n");

  fprintf(stderr, "num_scans: %d\n", cinfo->num_scans);
  fprintf(stderr, "arith_code: %d\n", cinfo->arith_code);
  fprintf(stderr, "optimize_coding: %d\n", cinfo->optimize_coding);
  fprintf(stderr, "CCIR601_sampling: %d\n", cinfo->CCIR601_sampling);
#if JPEG_LIB_VERSION >= 70
  fprintf(stderr, "do_fancy_downsampling: %d\n", cinfo->do_fancy_downsampling);
#endif
  fprintf(stderr, "smoothing_factor %d\n", cinfo->smoothing_factor);
  fprintf(stderr, "dct_method: %d\n", cinfo->dct_method);

  if (cinfo->master) {
    fprintf(stderr, "optimize_scans: %d\n", cinfo->master->optimize_scans);
    fprintf(stderr, "trellis_quant: %d\n", cinfo->master->trellis_quant);
    fprintf(stderr, "trellis_quant_dc: %d\n", cinfo->master->trellis_quant_dc);
    fprintf(stderr, "trellis_eob_opt: %d\n", cinfo->master->trellis_eob_opt);
    fprintf(stderr, "use_lambda_weight_tbl: %d\n", cinfo->master->use_lambda_weight_tbl);
    fprintf(stderr, "use_scans_in_trellis: %d\n", cinfo->master->use_scans_in_trellis);
    fprintf(stderr, "trellis_passes: %d\n", cinfo->master->trellis_passes);
    fprintf(stderr, "trellis_q_opt: %d\n", cinfo->master->trellis_q_opt);
    fprintf(stderr, "overshoot_deringing: %d\n", cinfo->master->overshoot_deringing);
    fprintf(stderr, "compress_profile: %d\n", cinfo->master->compress_profile);
    fprintf(stderr, "dc_scan_opt_mode: %d\n", cinfo->master->dc_scan_opt_mode);
    fprintf(stderr, "quant_tbl_master_idx: %d\n", cinfo->master->quant_tbl_master_idx);
    fprintf(stderr, "trellis_freq_split: %d\n", cinfo->master->trellis_freq_split);
    fprintf(stderr, "trellis_num_loops: %d\n", cinfo->master->trellis_num_loops);
    fprintf(stderr, "num_scans_luma: %d\n", cinfo->master->num_scans_luma);
    fprintf(stderr, "num_scans_luma_dc: %d\n", cinfo->master->num_scans_luma_dc);
    fprintf(stderr, "num_scans_chroma_dc: %d\n", cinfo->master->num_scans_chroma_dc);
    fprintf(stderr, "num_frequency_splits: %d\n", cinfo->master->num_frequency_splits);
    fprintf(stderr, "Al_max_luma: %d\n", cinfo->master->Al_max_luma);
    fprintf(stderr, "Al_max_chroma: %d\n", cinfo->master->Al_max_chroma);
    fprintf(stderr, "lambda_log_scale1: %f\n", cinfo->master->lambda_log_scale1);
    fprintf(stderr, "lambda_log_scale2: %f\n", cinfo->master->lambda_log_scale2);
    fprintf(stderr, "trellis_delta_dc_weight: %f\n", cinfo->master->trellis_delta_dc_weight);
  } else {
    fprintf(stderr, "cinfo->master was NULL\n");
  }
  fprintf(stderr, "\n");
}

*/
import "C"

import (
	"errors"
	"fmt"
	"image"
	"io"
	"unsafe"
)

// EncoderOptions specifies which settings to use during Compression.
type EncoderOptions struct {
	Quality         int
	OptimizeCoding  bool
	ProgressiveMode bool
	DCTMethod       DCTMethod
}

// Encode encodes src image and writes into w as JPEG format data.
func Encode(w io.Writer, src image.Image, opt *EncoderOptions) (err error) {
	// Recover panic
	defer func() {
		if r := recover(); r != nil {
			var ok bool
			err, ok = r.(error)
			if !ok {
				err = fmt.Errorf("JPEG error: %v", r)
			}
		}
	}()

	cinfo := C.new_compress()
	defer C.destroy_compress(cinfo)

	dstManager := makeDestinationManager(w, cinfo)
	defer releaseDestinationManager(dstManager)

	switch s := src.(type) {
	case *image.YCbCr:
		err = encodeYCbCr(cinfo, s, opt)
	case *image.Gray:
		err = encodeGray(cinfo, s, opt)
	case *image.NRGBA:
		err = encodeNRGBA(cinfo, s, opt)
	case *image.RGBA:
		err = encodeRGBA(cinfo, s, opt)
	default:
		return fmt.Errorf("unsupported image type: %T", src)
	}

	return
}

// encode image.YCbCr
func encodeYCbCr(cinfo *C.struct_jpeg_compress_struct, src *image.YCbCr, p *EncoderOptions) (err error) {
	// Set up compression parameters
	cinfo.image_width = C.JDIMENSION(src.Bounds().Dx())
	cinfo.image_height = C.JDIMENSION(src.Bounds().Dy())
	cinfo.input_components = 3
	cinfo.in_color_space = C.JCS_YCbCr

	C.jpeg_set_defaults(cinfo)
	setupEncoderOptions(cinfo, p)

	compInfo := (*[3]C.jpeg_component_info)(unsafe.Pointer(cinfo.comp_info))
	colorVDiv := 1
	switch src.SubsampleRatio {
	case image.YCbCrSubsampleRatio444:
		// 1x1,1x1,1x1
		compInfo[Y].h_samp_factor, compInfo[Y].v_samp_factor = 1, 1
		compInfo[Cb].h_samp_factor, compInfo[Cb].v_samp_factor = 1, 1
		compInfo[Cr].h_samp_factor, compInfo[Cr].v_samp_factor = 1, 1
	case image.YCbCrSubsampleRatio440:
		// 1x2,1x1,1x1
		compInfo[Y].h_samp_factor, compInfo[Y].v_samp_factor = 1, 2
		compInfo[Cb].h_samp_factor, compInfo[Cb].v_samp_factor = 1, 1
		compInfo[Cr].h_samp_factor, compInfo[Cr].v_samp_factor = 1, 1
		colorVDiv = 2
	case image.YCbCrSubsampleRatio422:
		// 2x1,1x1,1x1
		compInfo[Y].h_samp_factor, compInfo[Y].v_samp_factor = 2, 1
		compInfo[Cb].h_samp_factor, compInfo[Cb].v_samp_factor = 1, 1
		compInfo[Cr].h_samp_factor, compInfo[Cr].v_samp_factor = 1, 1
	case image.YCbCrSubsampleRatio420:
		// 2x2,1x1,1x1
		compInfo[Y].h_samp_factor, compInfo[Y].v_samp_factor = 2, 2
		compInfo[Cb].h_samp_factor, compInfo[Cb].v_samp_factor = 1, 1
		compInfo[Cr].h_samp_factor, compInfo[Cr].v_samp_factor = 1, 1
		colorVDiv = 2
	}

	// libjpeg raw data in is in planar format, which avoids unnecessary
	// planar->packed->planar conversions.
	cinfo.raw_data_in = C.TRUE

	// Start compression
	C.jpeg_start_compress(cinfo, C.TRUE)
	C.encode_ycbcr(
		cinfo,
		C.JSAMPROW(unsafe.Pointer(&src.Y[0])),
		C.JSAMPROW(unsafe.Pointer(&src.Cb[0])),
		C.JSAMPROW(unsafe.Pointer(&src.Cr[0])),
		C.int(src.YStride),
		C.int(src.CStride),
		C.int(colorVDiv),
	)

	C.print_cinfo(cinfo)

	C.jpeg_finish_compress(cinfo)
	return
}

// encode image.RGBA
func encodeRGBA(cinfo *C.struct_jpeg_compress_struct, src *image.RGBA, p *EncoderOptions) (err error) {
	// Set up compression parameters
	cinfo.image_width = C.JDIMENSION(src.Bounds().Dx())
	cinfo.image_height = C.JDIMENSION(src.Bounds().Dy())
	cinfo.input_components = 4
	cinfo.in_color_space = getJCS_EXT_RGBA()
	if cinfo.in_color_space == C.JCS_UNKNOWN {
		return errors.New("JCS_EXT_RGBA is not supported (probably built without libjpeg-turbo)")
	}

	C.jpeg_set_defaults(cinfo)
	setupEncoderOptions(cinfo, p)

	// Start compression
	C.jpeg_start_compress(cinfo, C.TRUE)
	C.encode_rgba(cinfo, C.JSAMPROW(unsafe.Pointer(&src.Pix[0])), C.int(src.Stride))

	C.print_cinfo(cinfo)

	C.jpeg_finish_compress(cinfo)
	return
}

// encode image.NRGBA
func encodeNRGBA(cinfo *C.struct_jpeg_compress_struct, src *image.NRGBA, p *EncoderOptions) (err error) {
	// Set up compression parameters
	cinfo.image_width = C.JDIMENSION(src.Bounds().Dx())
	cinfo.image_height = C.JDIMENSION(src.Bounds().Dy())
	cinfo.input_components = 4
	cinfo.in_color_space = getJCS_EXT_RGBA()
	if cinfo.in_color_space == C.JCS_UNKNOWN {
		return errors.New("JCS_EXT_RGBA is not supported (probably built without libjpeg-trubo)")
	}

	C.jpeg_set_defaults(cinfo)
	setupEncoderOptions(cinfo, p)

	// Start compression
	C.jpeg_start_compress(cinfo, C.TRUE)
	C.encode_rgba(cinfo, C.JSAMPROW(unsafe.Pointer(&src.Pix[0])), C.int(src.Stride))

	C.print_cinfo(cinfo)

	C.jpeg_finish_compress(cinfo)
	return
}

// encode image.Gray
func encodeGray(cinfo *C.struct_jpeg_compress_struct, src *image.Gray, p *EncoderOptions) (err error) {
	// Set up compression parameters
	cinfo.image_width = C.JDIMENSION(src.Bounds().Dx())
	cinfo.image_height = C.JDIMENSION(src.Bounds().Dy())
	cinfo.input_components = 1
	cinfo.in_color_space = C.JCS_GRAYSCALE

	C.jpeg_set_defaults(cinfo)
	setupEncoderOptions(cinfo, p)

	compInfo := (*C.jpeg_component_info)(unsafe.Pointer(cinfo.comp_info))
	compInfo.h_samp_factor, compInfo.v_samp_factor = 1, 1

	// libjpeg raw data in is in planar format, which avoids unnecessary
	// planar->packed->planar conversions.
	cinfo.raw_data_in = C.TRUE

	// Start compression
	C.jpeg_start_compress(cinfo, C.TRUE)
	C.encode_gray(cinfo, C.JSAMPROW(unsafe.Pointer(&src.Pix[0])), C.int(src.Stride))
	C.jpeg_finish_compress(cinfo)
	return
}

func setupEncoderOptions(cinfo *C.struct_jpeg_compress_struct, opt *EncoderOptions) {
	C.jpeg_set_quality(cinfo, C.int(opt.Quality), C.TRUE)
	if opt.OptimizeCoding {
		cinfo.optimize_coding = C.TRUE
	} else {
		cinfo.optimize_coding = C.FALSE
	}
	if opt.ProgressiveMode {
		C.jpeg_simple_progression(cinfo)
	}
	cinfo.dct_method = C.J_DCT_METHOD(opt.DCTMethod)

	C.print_cinfo(cinfo)
}
