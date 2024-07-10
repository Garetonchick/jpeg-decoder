#include <fft.h>

#include <fftw3.h>
#include <stdexcept>
#include <cmath>

class DctCalculator::Impl {
public:
    Impl(size_t width, std::vector<double> *input, std::vector<double> *output)
        : width_(width), input_(input), output_(output) {
        if (input_ == nullptr || output_ == nullptr || input_->size() != width_ * width_ ||
            output_->size() != width_ * width_) {
            throw std::invalid_argument("DctCalculator::Constructor: Invalid arguments");
        }
        plan_ = fftw_plan_r2r_2d(width_, width_, &(*input_)[0], &(*output_)[0],
                                 fftw_r2r_kind::FFTW_REDFT01, fftw_r2r_kind::FFTW_REDFT01,
                                 FFTW_MEASURE);
    }

    void Inverse() {
        for (size_t i = 0; i < width_ * width_; ++i) {
            double coef = 1.0;

            if (i < width_) {
                coef *= sqroot2_;
            }

            if (i % width_ == 0) {
                coef *= sqroot2_;
            }

            (*input_)[i] *= coef;
        }

        fftw_execute(plan_);

        for (size_t i = 0; i < width_ * width_; ++i) {
            (*output_)[i] /= 16.0;
        }
    }

    ~Impl() {
        fftw_destroy_plan(plan_);
    }

private:
    const double sqroot2_ = sqrt(2.0);

    size_t width_;
    std::vector<double> *input_;
    std::vector<double> *output_;
    fftw_plan plan_;
};

DctCalculator::DctCalculator(size_t width, std::vector<double> *input, std::vector<double> *output)
    : impl_(std::make_unique<DctCalculator::Impl>(width, input, output)) {
}

void DctCalculator::Inverse() {
    impl_->Inverse();
}

DctCalculator::~DctCalculator() = default;
